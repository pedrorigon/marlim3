# T082 — portão do Estágio 6

Corpus completo, 14 modelos. `MARLIM_DEFERRED_GATES=4`.

| gate | veredito | evidência |
|---|---|---|
| 1 — build limpo, sem aviso novo | **PASS** | 3 compilados, 0 erros, 162 avisos, 0 flags aumentadas |
| 2 — equivalência bit-a-bit (L2) | **PASS** | 14/14 modelos |
| 3 — suíte de regressão (L3) | **PASS** | 5 testes |
| 4 — desempenho | **DIFERIDO** | decisão do dono; **não verificado** |
| 5 — arquivos de referência intocados | **PASS** | `tests/comparison/` limpo |
| 6 — consumidores compilam | **PASS** | Num4Main, FA_Hidratos, FA_Hidratos_Servico, SisProdVap |

Resultado: **5 de 6 passaram, 1 diferido e não verificado.** O gate 4 reporta
DEFERRED, nunca PASS — o sumário não afirma cobertura que a corrida não tem.

## Modelos exigidos nominalmente pelo aceite

    2zones-2GLVs-2-Check-PA ......................... equivalent (1304 arquivos)
    extended-shutdown-combined-ESP-CGL-PIG-complete .. equivalent (102 arquivos)

Ambos com gas-lift, ambos bit-a-bit.

## Suplementar — calibração

Nove instrumentos, todos passando. Duas falhas foram encontradas e corrigidas
ao rodar (não ao ler) a bateria; estão descritas no commit do T081r. A segunda
é a mais séria: o `calibrate-all` pontuava uma sonda morta como `dead=0`, e só
saiu com erro por acidente do código de retorno do script filho.

## Suplementar — comparação estrutural L0

Primeira corrida: **22 falhas não declaradas**. São exatamente os 22 stubs de
delegação — o corpo saiu para `SisProdGasLift.cpp` e ficou uma chamada de uma
linha, então o L0 compara um corpo contra uma delegação e acusa DIFFERS em
todos.

Foram declaradas em `decomposed-functions.txt`, e a declaração é defensável
porque a supressão é estreita — verificado em `verify-structural.py`, não
assumido:

- `MISSING` conta como falha **incondicionalmente**, declarada ou não: uma
  função que desaparecer de todos os arquivos ainda falha;
- `STALE` dispara quando uma função declarada volta a ser idêntica, então uma
  entrada que sobrevive ao seu motivo se denuncia;
- cada corpo já fora provado token a token na sua própria tarefa (T073–T079), e
  qualquer mudança de comportamento cai no L2 e no gate 6.

Após a declaração: **0 falhas não declaradas, 71 declaradas.**

## Ressalva de ambiente

O pre-flight detectou que a **glibc mudou** desde a captura do baseline
(2.35-0ubuntu3.14 → 2.35-0ubuntu3.15). Isto já invalidou um baseline neste
projeto uma vez, silenciosamente, via `libmvec`.

O L2 passar 14/14 é, por si só, a evidência mais forte disponível: prova que
este ambiente, com o código atual, reproduz os artefatos do baseline bit a bit
— o que não seria possível se a atualização tivesse mexido nos resultados de
`libm` nestes caminhos. Ainda assim, `verify-baseline.sh` foi executado, porque
o harness manda verificar e não deduzir.

    BASELINE VALID -- an unmodified build still reproduces the capture
    equivalent models : 14/14

Um build do commit pristino `0f3b64f`, feito neste ambiente, reproduz a captura
nos 14 modelos. O baseline continua válido; a atualização de glibc não alcançou
estes caminhos. Registro completo em `verify-baseline.log`.
