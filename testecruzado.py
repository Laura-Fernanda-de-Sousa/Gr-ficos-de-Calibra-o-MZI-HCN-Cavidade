import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from scipy.optimize import curve_fit

# ==============================================================================
# 1. CARREGAMENTO E PREPARAÇÃO DOS DADOS
# ==============================================================================
dados = pd.read_csv("dados_mzi.csv")

# Pegamos apenas a ida (primeira metade)
meio = len(dados) // 2
transm = dados["Transm."].values[:meio]
mzi = dados["MZI"].values[:meio]
hcn = dados["HCN"].values[:meio]

# Identificação dos picos e vales 
# (Mantenha a proeminência do HCN em 0.05 a 0.08, conforme funcionou no seu diagnóstico)
picos_mzi, _ = find_peaks(mzi, prominence=0.05)
vales_mzi, _ = find_peaks(-mzi, prominence=0.05)
vales_hcn, _ = find_peaks(-hcn, prominence=0.05, distance=100)

extremos_pares_mzi = picos_mzi
extremos_impares_mzi = vales_mzi

# ==============================================================================
# 2. TABELA DO HCN (29 PICOS REAIS: R2 ao R0 + P1 ao P26)
# ==============================================================================
freqs_literatura = np.array([
    194615.93, 194533.11, 194449.11, 194277.55, 194190.04, 
    194101.35, 194011.50, 193920.50, 193828.34, 193735.02, 
    193640.54, 193544.91, 193448.12, 193350.19, 193251.10, 
    193150.87, 193049.50, 192946.98, 192843.33, 192738.53, 
    192632.60, 192525.54, 192417.35, 192308.03, 192197.58, 
    192086.02, 191973.34, 191859.54, 191744.63
])

# Trava de segurança para garantir tamanhos iguais
min_picos = min(len(vales_hcn), len(freqs_literatura))
vales_hcn_validos = vales_hcn[:min_picos]
freqs_lit_validas = freqs_literatura[:min_picos]

# Separação Pares (Calibração) e Ímpares (Teste)
hcn_para_calibrar = vales_hcn_validos[0::2] 
hcn_para_testar = vales_hcn_validos[1::2]   
freqs_calibracao_lit = freqs_lit_validas[0::2]
freqs_teste_lit = freqs_lit_validas[1::2]

# ==============================================================================
# 3. CALIBRAÇÃO ABSOLUTA E EIXO DE FREQUÊNCIAS
# ==============================================================================
indice_hcn_inicial = hcn_para_calibrar[0]
indice_hcn_final = hcn_para_calibrar[-1]
freq_inicial = freqs_calibracao_lit[0]
freq_final = freqs_calibracao_lit[-1]

# Identifica a direção da varredura (Alta -> Baixa frequência)
sinal_varredura = np.sign(freq_final - freq_inicial)

qtd_fsr_pares = np.sum((extremos_pares_mzi >= indice_hcn_inicial) & (extremos_pares_mzi <= indice_hcn_final))
FSR_pares = abs(freq_final - freq_inicial) / qtd_fsr_pares

qtd_fsr_impares = np.sum((extremos_impares_mzi >= indice_hcn_inicial) & (extremos_impares_mzi <= indice_hcn_final))
FSR_impares = abs(freq_final - freq_inicial) / qtd_fsr_impares

FSR_medio = (FSR_pares + FSR_impares) / 2
extremos_mzi = np.sort(np.concatenate((picos_mzi, vales_mzi)))

distancia_em_extremos = np.arange(len(extremos_mzi)) - np.searchsorted(extremos_mzi, indice_hcn_inicial)
freq_nos_extremos = freq_inicial + (distancia_em_extremos * sinal_varredura * (FSR_medio / 2))

indices_totais = np.arange(len(mzi))
eixo_freq_GHz = np.interp(indices_totais, extremos_mzi, freq_nos_extremos)

# ==============================================================================
# 4. VALIDAÇÃO CRUZADA (GRÁFICO PARA RESUMO)
# ==============================================================================
erros_validacao = []
for i, indice_teste in enumerate(hcn_para_testar):
    freq_medida = eixo_freq_GHz[indice_teste]
    freq_real = freqs_teste_lit[i]
    erros_validacao.append(freq_medida - freq_real)

plt.style.use('default')
plt.rcParams.update({'font.size': 12, 'axes.facecolor': 'white', 'figure.facecolor': 'white'})
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), gridspec_kw={'height_ratios': [3, 1]}, sharex=False)
fig.subplots_adjust(hspace=0.4) # Aumentei o espaço para caber o texto rotacionado

# Painel Superior
ax1.plot(eixo_freq_GHz, hcn, color='black', linewidth=1, label="Sinal HCN")
ax1.plot(eixo_freq_GHz[hcn_para_calibrar], hcn[hcn_para_calibrar], 'bo', markersize=8, label="Pares (Calibração)")
ax1.plot(eixo_freq_GHz[hcn_para_testar], hcn[hcn_para_testar], 'rX', markersize=8, label="Ímpares (Validação)")
ax1.set_ylabel("Transmissão HCN (V)")
ax1.set_title("Validação Cruzada da Calibração do Interferômetro")
ax1.legend(loc="best")
ax1.grid(True, linestyle='--', alpha=0.6)

# Painel Inferior (Erro)
ax2.axhline(0, color='black', linewidth=1) 
ax2.bar(range(len(erros_validacao)), erros_validacao, color='red', alpha=0.7, width=0.4)
ax2.set_xticks(range(len(erros_validacao)))
# Aqui as letras ficam menores e inclinadas para não grudarem:
ax2.set_xticklabels([f"Teste {i+1}" for i in range(len(erros_validacao))], fontsize=9, rotation=45, ha='right')
ax2.set_ylabel("Erro (GHz)")
ax2.set_xlabel("Picos de Teste")
ax2.set_title("Erro Absoluto (Medido vs Literatura)")
ax2.grid(True, axis='y', linestyle='--', alpha=0.6)

plt.tight_layout()
plt.savefig("validacao_hcn_resumo.png", dpi=300)

# ==============================================================================
# 5. AJUSTE NÃO LINEAR DA CAVIDADE (ROBUSTO)
# ==============================================================================
def lorentziana(f, f_c, gamma, A, B):
    return A - (B / ((f - f_c)**2 + (gamma/2)**2))

vales_transm, prop_transm = find_peaks(-transm, prominence=0.005) 

if len(vales_transm) == 0:
    print("\n[AVISO] Nenhuma ressonância clara encontrada no sinal de transmissão.")
else:
    indice_ressonancia = vales_transm[np.argmax(prop_transm['prominences'])]
    tamanho_janela = 2000 
    inicio_ressonancia = max(0, indice_ressonancia - tamanho_janela)
    fim_ressonancia = min(len(transm), indice_ressonancia + tamanho_janela)
    
    f_recorte = eixo_freq_GHz[inicio_ressonancia:fim_ressonancia]
    t_recorte = transm[inicio_ressonancia:fim_ressonancia]
    
    chute_fc = eixo_freq_GHz[indice_ressonancia]
    chute_gamma = 0.5 
    chute_A = np.mean(t_recorte[:100]) 
    chute_B = (chute_A - np.min(t_recorte)) * (chute_gamma/2)**2
    
    limite_inferior = [min(f_recorte), 0.0001, 0, 0]
    limite_superior = [max(f_recorte), 10.0, np.inf, np.inf]
    
    try:
        parametros_otimizados, _ = curve_fit(
            lorentziana, f_recorte, t_recorte, 
            p0=[chute_fc, chute_gamma, chute_A, chute_B], 
            bounds=(limite_inferior, limite_superior)
        )
        fc_ajuste, gamma_ajuste, A_ajuste, B_ajuste = parametros_otimizados
        fator_Q = fc_ajuste / gamma_ajuste
        
        plt.figure(figsize=(8, 5))
        plt.plot(f_recorte, t_recorte, 'k.', markersize=2, label="Dados Experimentais")
        plt.plot(f_recorte, lorentziana(f_recorte, *parametros_otimizados), 'r-', linewidth=2, label=f"Ajuste (Q = {fator_Q:.1e})")
        plt.xlabel("Frequência (GHz)")
        plt.ylabel("Transmissão (V)")
        plt.title("Ajuste Não Linear da Ressonância Óptica")
        plt.legend()
        plt.grid(True, linestyle='--', alpha=0.6)
        
    except RuntimeError:
        print("\n[ERRO] O curve_fit não conseguiu convergir para um resultado físico.")

# ==============================================================================
# 6. PAINEL DE RESULTADOS EM IMAGEM (A "IMAGEM BONITINHA")
# ==============================================================================
plt.figure(figsize=(5, 3), facecolor='white')
plt.axis('off')

# Cálculo do erro médio para mostrar no quadro
erro_medio = np.mean(np.abs(erros_validacao))

texto = (
    "--- TESTE DE SIMETRIA DO MZI ---\n"
    f"FSR (Pares): {FSR_pares:.4f} GHz\n"
    f"FSR (Ímpares): {FSR_impares:.4f} GHz\n"
    f"Diferença: {abs(FSR_pares - FSR_impares):.6f} GHz\n\n"
    "--- VALIDAÇÃO CRUZADA ---\n"
    f"Erro Médio: {erro_medio:.4f} GHz"
)

plt.text(0.5, 0.5, texto, fontsize=11, family='monospace', 
         ha='center', va='center', 
         bbox=dict(boxstyle='round,pad=1', facecolor='#f8f9fa', edgecolor='#333333', linewidth=1.5))

plt.tight_layout()
plt.savefig("simetria_mzi_resumo.png", dpi=300, bbox_inches='tight')

# Abre todas as imagens geradas na tela!
plt.show()