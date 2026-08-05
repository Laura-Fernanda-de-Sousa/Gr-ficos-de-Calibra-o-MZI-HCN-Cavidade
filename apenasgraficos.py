#Importar bibliotecas importantes
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.signal import find_peaks

#Carregar os dados de forma que utilizemos apenas a ida
dados = pd.read_csv("dados_mzi.csv")
meio = len(dados) // 2  #varredura completa dividida por 2
transm = dados["Transm."].values[:meio]
mzi = dados["MZI"].values[:meio]
hcn = dados["HCN"].values[:meio]

#Identificar picos e vales
picos_mzi, _ = find_peaks(mzi, prominence=0.05)
vales_mzi, _ = find_peaks(-mzi, prominence=0.05)
extremos_mzi = np.sort(np.concatenate((picos_mzi, vales_mzi)))
vales_hcn, _ = find_peaks(-hcn, prominence=0.1, distance=100)  #Encontrar vales do HCN

#Calibrar MZI em relacao ao HCN
freq_hcn_1_GHz = 193100.0  #Ajustado para ~193100 GHz para uma leitura mais refinada
indice_hcn_1 = vales_hcn[0]
freq_hcn_2_GHz = 193500.0 
indice_hcn_2 = vales_hcn[-1]
extremos_entre_referencias = np.sum((extremos_mzi >= indice_hcn_1) & (extremos_mzi <= indice_hcn_2))
delta_freq = abs(freq_hcn_2_GHz - freq_hcn_1_GHz)
FSR_GHz = (delta_freq / extremos_entre_referencias) * 2
print(f"Calibração Concluída no sinal de ida!")
print(f"FSR Calculado: {FSR_GHz:.4f} GHz")

#Estabelecer o eixo de frequencia
distancia_em_extremos = np.arange(len(extremos_mzi)) - np.searchsorted(extremos_mzi, indice_hcn_1)
freq_nos_extremos = freq_hcn_1_GHz + (distancia_em_extremos * (FSR_GHz / 2))
indices_totais = np.arange(len(mzi))
eixo_freq_GHz = np.interp(indices_totais, extremos_mzi, freq_nos_extremos)

#Plotar os graficos de forma independente
plt.style.use('default') #Fundo branco e escrita preta
plt.rcParams.update({
    'axes.facecolor': 'white',
    'figure.facecolor': 'white',
    'text.color': 'black',
    'axes.labelcolor': 'black',
    'xtick.color': 'black',
    'ytick.color': 'black',
    'grid.color': 'gray',
    'grid.alpha': 0.5
})

#GRAFICO 1: TRANSMISSAO
plt.figure(1, figsize=(8, 5))
plt.plot(eixo_freq_GHz, transm, color='blue', linewidth=1.5, label="Transmissão")
plt.xlabel("Frequência (GHz)")
plt.ylabel("Intensidade (V)")
plt.title("Cavidade Fotônica")
plt.legend(loc="best")
plt.grid(True, linestyle='--')
plt.tight_layout()

#GRAFICO 2: MZI
plt.figure(2, figsize=(8, 5))
plt.plot(eixo_freq_GHz, mzi, color='blue', linewidth=1.5, label=f"FSR = {FSR_GHz:.2f} GHz")
plt.xlabel("Frequência (GHz)")
plt.ylabel("Intensidade (V)")
plt.title("Interferômetro Mach-Zehnder")
plt.legend(loc="best")
plt.grid(True, linestyle='--')
plt.tight_layout()

#GRAFICO 3: HCN
plt.figure(3, figsize=(8, 5))
plt.plot(eixo_freq_GHz, hcn, color='blue', linewidth=1.5, label="Espectro HCN")
plt.plot(eixo_freq_GHz[indice_hcn_1], hcn[indice_hcn_1], 'or', markersize=6, label="Ref. 1")
plt.plot(eixo_freq_GHz[indice_hcn_2], hcn[indice_hcn_2], 'ob', markersize=6, label="Ref. 2")
plt.xlabel("Frequência (GHz)")
plt.ylabel("Intensidade (V)")
plt.title("Gás de Referência")
plt.legend(loc="best")
plt.grid(True, linestyle='--')
plt.tight_layout()

#GRAFICO 4: CURVA DE CALIBRACAO--
plt.figure(4, figsize=(8, 5))
plt.plot(indices_totais, eixo_freq_GHz, color='blue', linewidth=1.5, label="Interpolação Linear")
plt.plot(extremos_mzi, freq_nos_extremos, 'o', color='red', markersize=4, label="Pontos MZI (FSR/2)")
plt.xlabel("Índice da Amostra (Tempo)")
plt.ylabel("Frequência (GHz)")
plt.title("Curva de Calibração: Frequência vs Tempo")
plt.legend(loc="best")
plt.grid(True, linestyle='--')
plt.tight_layout()

#Abrir as 4 janelas ao mesmo tempo
plt.show()