#!/usr/bin/env bash
# Calibrate the time-step policy expansion in transient-solve-move.py (T109a).
#
# T109a replaced three conditional calls in the transient step with a registry
# of policies. The token proof of the step now expands every call to the
# registry from the policy definitions in the module -- which is only worth
# something if a wrong definition makes the proof fail. This script shows that
# it does: each case below corrupts one thing the expansion reads, in a COPY of
# the module, and the check must reject it.
#
# The corruptions are the ways a registry goes wrong without anything else
# noticing: a changed condition, a changed action, a policy moved to another
# hook, a policy dropped from the list or listed twice, the call's context bound
# in the wrong order, a delegation to the wrong policy, and a policy the
# expansion cannot read at all. Every one of them compiles.
#
# Three controls must pass: the module untouched, for each of the two functions
# that call the registry, and a test policy that only delegates to an existing
# one, in that one's place in the list. The last is T109a's own acceptance case,
# and it passing shows the expansion resolves delegation instead of rejecting
# everything it did not write.
#
# Each case asserts the corruption was injected before running the checker; a
# replacement that matched nothing would test nothing.
#
# The module is read from the commit that introduced the registry, so renames
# made later cannot erode this instrument, the way they eroded others. Override
# with MARLIM_T109A_COMMIT, or point MARLIM_T109A_MODULE at a file.
#
# Usage: calibrate-time-step-policies.sh
# Exit code: 0 when every control passes and every corruption is caught.

set -uo pipefail
export LC_ALL=C

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
project_root="$(cd "$script_dir/.." && pwd)"

PRE_T127="${MARLIM_PRE_T127:-baa24f0^}"
T109A_COMMIT="${MARLIM_T109A_COMMIT:-0828baa}"

red=$'\033[0;31m'; green=$'\033[0;32m'; bold=$'\033[1m'; reset=$'\033[0m'

work="$(mktemp -d -t marlim3-tsp-cal-XXXXXX)"
trap 'rm -rf "$work"' EXIT

git -C "$project_root" show "$PRE_T127:src/core/SisProd.cpp" > "$work/baseline.cpp" || exit 2
if [[ -n "${MARLIM_T109A_MODULE:-}" ]]; then
    cp "$MARLIM_T109A_MODULE" "$work/module.cpp" || exit 2
else
    git -C "$project_root" show "$T109A_COMMIT:src/core/SisProdTransient.cpp" \
        > "$work/module.cpp" || exit 2
fi

failures=0
cases=0

# $1 label   $2 function to check (advanceCouplingIteration | SolveTrans)
# $3 python replacements, as a list of (old, new) pairs; empty for the module
#    untouched   $4 "control" when the replacements must still pass the check
attempt() {
    local label="$1" function="$2" replacements="$3" expect="${4:-}"
    [[ -z "$replacements" ]] && expect="control"
    cases=$((cases + 1))
    cp "$work/module.cpp" "$work/case.cpp"
    if [[ -n "$replacements" ]]; then
        if ! python3 - "$work/case.cpp" <<EOF
import sys
path = sys.argv[1]
text = open(path, encoding="utf-8").read()
for old, new in $replacements:
    if text.count(old) != 1:
        sys.exit(1)
    text = text.replace(old, new)
open(path, "w", encoding="utf-8").write(text)
EOF
        then
            printf '%s%-58s CORRUPTION NOT INJECTED%s\n' "$red" "$label" "$reset"
            printf '   a replacement matched zero or several times; this case tested nothing\n'
            failures=$((failures + 1))
            return
        fi
    fi

    local out status
    out="$(python3 "$script_dir/transient-solve-move.py" check "$function" \
             "$work/baseline.cpp" "$work/case.cpp" 2>&1)"
    status=$?

    if [[ "$expect" == "control" ]]; then
        if (( status == 0 )); then
            printf '%s%-58s PASSES (control)%s\n' "$green" "$label" "$reset"
        else
            printf '%s%-58s CONTROL FAILED%s\n' "$red" "$label" "$reset"
            printf '%s\n' "$out" | tail -4 | sed 's/^/   /'
            failures=$((failures + 1))
        fi
        return
    fi

    if (( status != 0 )) && ! grep -q '^OK ' <<< "$out"; then
        printf '%s%-58s caught%s   %s\n' "$green" "$label" "$reset" \
               "$(grep -m1 -oE '^DIFFERS.*|[A-Za-z]*Error: .{0,60}' <<< "$out" || true)"
    else
        printf '%s%-58s NOT CAUGHT%s\n' "$red" "$label" "$reset"
        failures=$((failures + 1))
    fi
}

printf '%s=== time-step policy expansion calibration ===%s\n\n' "$bold" "$reset"

attempt "control: module untouched (coupling iteration)" advanceCouplingIteration ""
attempt "control: module untouched (step)" SolveTrans ""

# T109a's acceptance case: a test policy that only delegates, in the place of
# the policy it delegates to. The step does what it did, so the proof passes.
TEST_POLICY='
struct TestDelegatingTimeStepPolicy {
    static constexpr TimeStepHook hook = DampMaximumTimeStepPolicy::hook;
    static bool applies(const TimeStepPolicyContext &context) {
        return DampMaximumTimeStepPolicy::applies(context);
    }
    static void apply(const TimeStepPolicyContext &context) {
        DampMaximumTimeStepPolicy::apply(context);
    }
};
'
attempt "control: delegating test policy in the delegate's place" advanceCouplingIteration \
    "[('/// The policies of a hook, applied', '''$TEST_POLICY
/// The policies of a hook, applied'''),
      ('    DampMaximumTimeStepPolicy,\n    PressureRateOfChangePolicy>;',
       '    TestDelegatingTimeStepPolicy,\n    PressureRateOfChangePolicy>;')]" control

attempt "condition: couplingIteration == 0 -> == 1" advanceCouplingIteration \
    "[('return context.couplingIteration == 0 &&', 'return context.couplingIteration == 1 &&')]"
attempt "condition: controleDTvalv == 1 -> == 0" advanceCouplingIteration \
    "[('context.step.input.controleDTvalv == 1;', 'context.step.input.controleDTvalv == 0;')]"
attempt "condition: 1 * fullModel -> fullModel" advanceCouplingIteration \
    "[('context.couplingIteration == 1 * context.step.fullModel;', 'context.couplingIteration == context.step.fullModel;')]"
attempt "action: evaluatePressureRateOfChange(..., 0, 0, ...) -> 0, 1" SolveTrans \
    "[('evaluatePressureRateOfChange(context.step, 0, 0, context.explicitScheme);', 'evaluatePressureRateOfChange(context.step, 0, 1, context.explicitScheme);')]"
attempt "action: dampMaximumTimeStep -> restrictTimeStepByValve" advanceCouplingIteration \
    "[('        dampMaximumTimeStep(context.step);', '        restrictTimeStepByValve(context.step);')]"
attempt "hook: damping moved to the start of the iteration" advanceCouplingIteration \
    "[('static constexpr TimeStepHook hook = TimeStepHook::AfterPigUpdate;', 'static constexpr TimeStepHook hook = TimeStepHook::CouplingIterationStart;')]"
attempt "hook: pressure rate moved into the coupling loop" SolveTrans \
    "[('static constexpr TimeStepHook hook = TimeStepHook::AfterValveOpenings;', 'static constexpr TimeStepHook hook = TimeStepHook::AfterPigUpdate;')]"
attempt "list: pressure-rate policy dropped" SolveTrans \
    "[('    DampMaximumTimeStepPolicy,\n    PressureRateOfChangePolicy>;', '    DampMaximumTimeStepPolicy>;')]"
attempt "list: damping listed twice" advanceCouplingIteration \
    "[('    DampMaximumTimeStepPolicy,\n', '    DampMaximumTimeStepPolicy,\n    DampMaximumTimeStepPolicy,\n')]"
attempt "list: delegating test policy added beside its delegate" advanceCouplingIteration \
    "[('/// The policies of a hook, applied', '''$TEST_POLICY
/// The policies of a hook, applied'''),
      ('    DampMaximumTimeStepPolicy,\n', '    DampMaximumTimeStepPolicy,\n    TestDelegatingTimeStepPolicy,\n')]"
attempt "context: couplingIteration and explicitScheme swapped" advanceCouplingIteration \
    "[('TimeStepPolicies::apply<TimeStepHook::AfterPigUpdate>({state.step, kontaAcop, vExpli});', 'TimeStepPolicies::apply<TimeStepHook::AfterPigUpdate>({state.step, vExpli, kontaAcop});')]"
attempt "delegation: test policy's condition from the wrong policy" advanceCouplingIteration \
    "[('/// The policies of a hook, applied', '''$TEST_POLICY
/// The policies of a hook, applied'''.replace('return DampMaximumTimeStepPolicy::applies', 'return RestrictTimeStepByValvePolicy::applies')),
      ('    DampMaximumTimeStepPolicy,\n    PressureRateOfChangePolicy>;',
       '    TestDelegatingTimeStepPolicy,\n    PressureRateOfChangePolicy>;')]"
attempt "form: a policy the expansion cannot read" advanceCouplingIteration \
    "[('    static void apply(const TimeStepPolicyContext &context) {\n        dampMaximumTimeStep(context.step);', '    static void apply(const TimeStepPolicyContext &) {\n        dampMaximumTimeStep(state.step);')]"

printf '\n'
if (( failures > 0 )); then
    printf '%sCALIBRATION FAILED -- %d of %d cases%s\n' "$red" "$failures" "$cases" "$reset" >&2
    exit 1
fi
printf '%sCALIBRATION PASSED -- %d cases: every control passes, every corruption caught%s\n' \
    "$green" "$cases" "$reset"
