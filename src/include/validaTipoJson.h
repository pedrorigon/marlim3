/*
 * validaJson.h
 *
 *  Created on: 20 de mar. de 2026
 *      Author: Eduardo
 */

#ifndef VALIDATIPOJSON_H_
#define VALIDATIPOJSON_H_

#include <string>
#include <vector>
using namespace std;
///////////////////////////////////////
#include "rapidjson/document.h"
#include "rapidjson/stringbuffer.h"
#include "rapidjson/writer.h"
#include "rapidjson/schema.h"
#include "rapidjson/filereadstream.h"
#include "rapidjson/filewritestream.h"
#include "rapidjson/error/pt_BR.h"
#include "JSON_entrada.h"
#include "Log.h"


#ifdef linux
  // barra para linux
  #define BARRA "/"
#elif LINUX
  // barra para linux
  #define BARRA "/"
#elif Linux
  // barra para linux
  #define BARRA "/"
#elif UNIX
  // barra para linux
  #define BARRA "/"
#elif defined WIN32
  // barra para windows
  #define BARRA "\\"
#elif defined Win32
  // barra para windows
  #define BARRA "\\"
#elif defined win32
  // barra para windows
  #define BARRA "\\"
#else
  // barra para outros
  #define BARRA "/\\"
#endif

using namespace rapidjson;

// relatorio dos arquivos de dados de saida da simulacao
extern ofstream arqRelatorioPerfis;

// criar objeto para log da aplicacao
extern Logger logger;

// obter string do path dos arquivos de entrada para o simulador
extern string pathArqExtEntrada;

// obter string do prefixo dos arquivos de saida para POCO_INJETOR
extern string pathPrefixoArqSaida;

extern string versao;


class validadorTipo{
	public:


	validadorTipo(){}

	void validaGeral(Document& jsonDoc);


    // Ma todos para validacao do arquivo de entrada
    void valida_configuracao_inicial(Value& configuracao_inicial_json, std::vector<std::string>& erros, bool& sucesso);
    void valida_condcont_pocinjec(Value& condcont_json, std::vector<std::string>& erros, bool& sucesso);
    void valida_tempo(Value& tempo_json, std::vector<std::string>& erros, bool& sucesso);
    void valida_materiais(Value& material_json, std::vector<std::string>& erros, bool& sucesso);
    void valida_corte(Value& corte_json, std::vector<std::string>& erros, bool& sucesso);
    void valida_tabela(Value& tabela_json, std::vector<std::string>& erros, bool& sucesso);
    void valida_parafina(Value& parafina_json, std::vector<std::string>& erros, bool& sucesso);
    void valida_fluidos_producao(Value& fluidos_producao_json, std::vector<std::string>& erros, bool& sucesso);
    void valida_fluido_complementar(Value& fluido_complementar_json, std::vector<std::string>& erros, bool& sucesso);
    void valida_fluido_gas(Value& fluido_gas_json, std::vector<std::string>& erros, bool& sucesso);

    void valida_unidades_producao(Value& duto_producao_json, std::vector<std::string>& erros, bool& sucesso);
    void valida_unidades_servico(Value& duto_servico_json, std::vector<std::string>& erros, bool& sucesso);
    void valida_hidrato(Value& hidrato_json, std::vector<std::string>& erros, bool& sucesso);
    void valida_ipr(Value& ipr_json, std::vector<std::string>& erros, bool& sucesso);
    void valida_separador(Value& separador_json, std::vector<std::string>& erros, bool& sucesso);
    void valida_correcao(Value& correcao_json, std::vector<std::string>& erros, bool& sucesso);
    void valida_chokeSup(Value& chokeSup_json, std::vector<std::string>& erros, bool& sucesso);
    void valida_chokeInj(Value& chokeInj_json, std::vector<std::string>& erros, bool& sucesso);
    void valida_master1(Value& master1_json, std::vector<std::string>& erros, bool& sucesso);
    void valida_master2(Value& master2_json, std::vector<std::string>& erros, bool& sucesso);
    void valida_valv(Value& valvula_json, std::vector<std::string>& erros, bool& sucesso);
    void valida_gasInj(Value& gasInj_json, std::vector<std::string>& erros, bool& sucesso);
    void valida_fonte_gas(Value& fonte_gas_json, std::vector<std::string>& erros, bool& sucesso);
    void valida_fonte_liquido(Value& fonte_liquido_json, std::vector<std::string>& erros, bool& sucesso);
    void valida_fonte_PoroRadial(Value& fonte_poroRadial_json, std::vector<std::string>& erros, bool& sucesso);
    void valida_fonte_Poro2D(Value& fonte_poro2D_json, std::vector<std::string>& erros, bool& sucesso);
    void valida_fonte_massa(Value& fonte_massa_json, std::vector<std::string>& erros, bool& sucesso);
    void valida_furo(Value& fontePressao_json, std::vector<std::string>& erros, bool& sucesso);
    void valida_fontechk(Value& fontechk_json, std::vector<std::string>& erros, bool& sucesso);
    void valida_pig(Value& pig_json, std::vector<std::string>& erros, bool& sucesso);
    void valida_bcs(Value& bcs_json, std::vector<std::string>& erros, bool& sucesso);
    void valida_multibcs(Value& multibcs_json, std::vector<std::string>& erros, bool& sucesso);
    void valida_bomba_volumetrica(Value& bomba_volumetrica_json, std::vector<std::string>& erros, bool& sucesso);
    void valida_delta_pressao(Value& delta_pressao_json, std::vector<std::string>& erros, bool& sucesso);
    void valida_fonteCalor(Value& fonteCalor_json, std::vector<std::string>& erros, bool& sucesso);
    void valida_fonte_gaslift(Value& fonte_gaslift_json, std::vector<std::string>& erros, bool& sucesso);
    void valida_intermitencia(Value& intermitencia_json, std::vector<std::string>& erros, bool& sucesso);

    void valida_perfil_producao(Value& perfis_producao_json, std::vector<std::string>& erros, bool& sucesso);
    void valida_perfil_servico(Value& perfis_servico_json, std::vector<std::string>& erros, bool& sucesso);

    void valida_tendencia_producao(Value& tendencia_producao_json, std::vector<std::string>& erros, bool& sucesso);
    void valida_tendencia_servico(Value& tendencia_servico_json, std::vector<std::string>& erros, bool& sucesso);




};



#endif /* VALIDATIPOJSON_H_ */
