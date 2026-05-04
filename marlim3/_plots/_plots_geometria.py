import math
import plotly.graph_objects as go


def _get_duto_length(duto):
    """Calculate duct length from discretization or dxCelula."""
    agrupamento = duto.get("agrupamento", True)
    if agrupamento:
        disc = duto.get("discretizacao", [])
        if not isinstance(disc, list) or not disc:
            return 0
        return sum(
            (item.get("nCelulas", 0) or 0) * (item.get("comprimento", 0) or 0)
            for item in disc if isinstance(item, dict)
        )
    else:
        dx_celula = duto.get("dxCelula", [])
        if not isinstance(dx_celula, list) or not dx_celula:
            return 0
        return sum(float(x) for x in dx_celula)


def _plotar_geometria(tramo):

    # Configurações iniciais
    x_coords_prod = [0]  # Coordenada X inicial para dutos de produção
    y_coords_prod = [0]  # Coordenada Y inicial para dutos de produção
    tooltips_prod = []  # Tooltip para dutos de produção

    if tramo.dutosServico:
        x_coords_serv = []  # Coordenada X inicial para dutos de serviço
        y_coords_serv = []  # Coordenada Y inicial para dutos de serviço
        tooltips_serv = []  # Tooltip para dutos de serviço

    # Determine if modoXY is active
    config_inicial = getattr(tramo, "configuracaoInicial", None) or {}
    if isinstance(config_inicial, dict):
        modo_xy = config_inicial.get("modoXY", False)
    else:
        modo_xy = getattr(config_inicial, "modoXY", False)

    # Processando dutos de produção
    for idx, duto in enumerate(tramo.dutosProducao):
        if modo_xy and "xCoor" in duto and "yCoor" in duto:
            # In XY mode, xCoor/yCoor represent absolute endpoint coordinates.
            # For production ducts starting from origin, we can use them directly.
            x_coords_prod.append(float(duto["xCoor"]))
            y_coords_prod.append(float(duto["yCoor"]))
        elif "angulo" in duto:
            # Angle mode: compute displacement from angle and length
            comprimento = _get_duto_length(duto)
            ang = float(duto["angulo"])
            dx = comprimento * math.cos(ang)
            dy = comprimento * math.sin(ang)
            x_coords_prod.append(x_coords_prod[-1] + dx)
            y_coords_prod.append(y_coords_prod[-1] + dy)
        else:
            continue
        tooltips_prod.append(
            f"ID: {duto.get('id', '?')}<br>"
            f"Rótulo: {duto.get('rotulo', 'N/A')}<br>"
            f"idCorte: {duto.get('idCorte', 'N/A')}<br>"
            f"Acoplamento Térmico: {duto.get('acoplamentoTermico', 'N/A')}"
        )

    # Coordenadas iniciais da plataforma (último ponto da produção)
    platform_x = x_coords_prod[-1]
    platform_y = y_coords_prod[-1]

    if tramo.dutosServico:
        # Coordenadas iniciais para os dutos de serviço
        # (começam na plataforma)
        x_coords_serv.append(platform_x)
        y_coords_serv.append(platform_y)

    # Processando dutos de serviço
    if tramo.dutosServico:
        for idx, duto in enumerate(tramo.dutosServico):
            if modo_xy and "xCoor" in duto and "yCoor" in duto:
                # In XY mode, compute delta from previous duct's xCoor/yCoor
                # (first duct's reference is origin 0,0)
                if idx == 0:
                    prev_x, prev_y = 0.0, 0.0
                else:
                    prev_duto = tramo.dutosServico[idx - 1]
                    prev_x = float(prev_duto.get("xCoor", 0))
                    prev_y = float(prev_duto.get("yCoor", 0))
                dx = float(duto["xCoor"]) - prev_x
                dy = float(duto["yCoor"]) - prev_y
                x_coords_serv.append(x_coords_serv[-1] + dx)
                y_coords_serv.append(y_coords_serv[-1] + dy)
            elif "angulo" in duto:
                # Angle mode: invert X direction for service ducts
                comprimento = _get_duto_length(duto)
                ang = float(duto["angulo"])
                dx = comprimento * math.cos(ang)
                dy = comprimento * math.sin(ang)
                x_coords_serv.append(x_coords_serv[-1] - dx)
                y_coords_serv.append(y_coords_serv[-1] + dy)
            else:
                continue
            tooltips_serv.append(
                f"ID: {duto.get('id', '?')}<br>"
                f"Rótulo: {duto.get('rotulo', 'N/A')}<br>"
                f"idCorte: {duto.get('idCorte', 'N/A')}<br>"
                f"Acoplamento Térmico: {duto.get('acoplamentoTermico', 'N/A')}"
            )

    # Criando o gráfico interativo
    fig = go.Figure()

    # Adicionando os dutos de produção
    fig.add_trace(go.Scatter(
        x=x_coords_prod,
        y=y_coords_prod,
        mode="lines+markers",
        name="Production",
        hovertext=tooltips_prod,
        hoverinfo="text",
        line=dict(color="#39C0E0", width=3),
        marker=dict(size=8, color="#39C0E0")
    ))

    # Adicionando os dutos de serviço
    if tramo.dutosServico:
        fig.add_trace(go.Scatter(
            x=x_coords_serv,
            y=y_coords_serv,
            mode="lines+markers",
            name="Service",
            hovertext=tooltips_serv,
            hoverinfo="text",
            line=dict(color="#FFA933", width=3),
            marker=dict(size=8, color="#FFA933")
        ))

    # Configurações de layout
    fig.update_layout(
        xaxis_title="x (m)",
        yaxis_title="y (m)",
        showlegend=False,
        xaxis=dict(scaleanchor="y"),
        yaxis=dict(scaleanchor="x")
    )

    # Mostrando o gráfico
    fig.show()
