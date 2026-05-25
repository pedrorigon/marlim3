import math
import plotly.graph_objects as go


def _g(d, *keys, default=None):
    """Get value from dict trying multiple key names (PT/EN compatibility)."""
    for k in keys:
        if k in d:
            return d[k]
    return default


def _get_duto_length(duto):
    """Calculate duct length from discretization or dxCelula."""
    agrupamento = _g(duto, "agrupamento", "grouping", default=True)
    if agrupamento:
        disc = _g(duto, "discretizacao", "discretization", default=[])
        if not isinstance(disc, list) or not disc:
            return 0
        return sum(
            (_g(item, "nCelulas", "numCells", default=0) or 0) *
            (_g(item, "comprimento", "length", default=0) or 0)
            for item in disc if isinstance(item, dict)
        )
    else:
        dx_celula = _g(duto, "dxCelula", "cellDx", default=[])
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
        modo_xy = _g(config_inicial, "modoXY", "xyMode", default=False)
        segue_escoamento = _g(config_inicial, "sentidoGeometriaSegueEscoamento",
                              "geometryFollowsFlow", default=True)
    else:
        modo_xy = getattr(config_inicial, "modoXY", False)
        segue_escoamento = getattr(config_inicial, "sentidoGeometriaSegueEscoamento", True)

    # When sentidoGeometriaSegueEscoamento=false:
    # - Duct 0 is at the platform, last duct is at the reservoir
    # - Angles are defined relative to flow direction (reservoir→platform)
    # - For geometry plotting, we reverse ducts to traverse reservoir→platform
    # - The -dx convention in service already mirrors the x-direction
    if segue_escoamento:
        dutos_prod_plot = tramo.dutosProducao
        angle_offset_prod = 0.0
        angle_offset_serv = 0.0
    else:
        dutos_prod_plot = list(reversed(tramo.dutosProducao))
        angle_offset_prod = 0.0
        angle_offset_serv = 0.0

    # Processando dutos de produção
    for idx, duto in enumerate(dutos_prod_plot):
        if modo_xy and (_g(duto, "xCoor", "xCoord") is not None) and (_g(duto, "yCoor", "yCoord") is not None):
            x_coords_prod.append(float(_g(duto, "xCoor", "xCoord")))
            y_coords_prod.append(float(_g(duto, "yCoor", "yCoord")))
        elif _g(duto, "angulo", "angle") is not None:
            # Angle mode: compute displacement from angle and length
            comprimento = _get_duto_length(duto)
            ang = float(_g(duto, "angulo", "angle")) + angle_offset_prod
            dx = comprimento * math.cos(ang)
            dy = comprimento * math.sin(ang)
            x_coords_prod.append(x_coords_prod[-1] + dx)
            y_coords_prod.append(y_coords_prod[-1] + dy)
        else:
            continue
        tooltips_prod.append(
            f"ID: {_g(duto, 'id', default='?')}<br>"
            f"Rótulo: {_g(duto, 'rotulo', 'label', default='N/A')}<br>"
            f"idCorte: {_g(duto, 'idCorte', 'crossSectionId', default='N/A')}<br>"
            f"Acoplamento Térmico: {_g(duto, 'acoplamentoTermico', 'thermalCoupling', default='N/A')}"
        )

    # Coordenadas iniciais da plataforma
    # In both cases, after traversal production ends at the platform
    platform_x = x_coords_prod[-1]
    platform_y = y_coords_prod[-1]

    if tramo.dutosServico:
        # Coordenadas iniciais para os dutos de serviço
        # (começam na plataforma)
        x_coords_serv.append(platform_x)
        y_coords_serv.append(platform_y)

    # Processando dutos de serviço
    # Service flows platform→reservoir; when segue_escoamento=false, duct 0
    # is already at the platform (our starting point), so no reversal needed.
    if tramo.dutosServico:
        dutos_serv_plot = tramo.dutosServico
        for idx, duto in enumerate(dutos_serv_plot):
            if modo_xy and (_g(duto, "xCoor", "xCoord") is not None) and (_g(duto, "yCoor", "yCoord") is not None):
                if idx == 0:
                    prev_x, prev_y = 0.0, 0.0
                else:
                    prev_duto = dutos_serv_plot[idx - 1]
                    prev_x = float(_g(prev_duto, "xCoor", "xCoord", default=0))
                    prev_y = float(_g(prev_duto, "yCoor", "yCoord", default=0))
                dx = float(_g(duto, "xCoor", "xCoord")) - prev_x
                dy = float(_g(duto, "yCoor", "yCoord")) - prev_y
                x_coords_serv.append(x_coords_serv[-1] + dx)
                y_coords_serv.append(y_coords_serv[-1] + dy)
            elif _g(duto, "angulo", "angle") is not None:
                comprimento = _get_duto_length(duto)
                ang = float(_g(duto, "angulo", "angle")) + angle_offset_serv
                dx = comprimento * math.cos(ang)
                dy = comprimento * math.sin(ang)
                x_coords_serv.append(x_coords_serv[-1] - dx)
                y_coords_serv.append(y_coords_serv[-1] + dy)
            else:
                continue
            tooltips_serv.append(
                f"ID: {_g(duto, 'id', default='?')}<br>"
                f"Rótulo: {_g(duto, 'rotulo', 'label', default='N/A')}<br>"
                f"idCorte: {_g(duto, 'idCorte', 'crossSectionId', default='N/A')}<br>"
                f"Acoplamento Térmico: {_g(duto, 'acoplamentoTermico', 'thermalCoupling', default='N/A')}"
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
