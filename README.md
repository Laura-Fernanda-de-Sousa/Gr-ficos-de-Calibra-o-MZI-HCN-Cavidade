# Graficos-de-Calibracao-MZI-HCN-Cavidade:

  O primeiro código, nomeado de "apenasgraficos.py", contém os comandos necessários para gerar a Calibração do MZI, os dados brutos do MZI,   dados do HCN para tal calibração e da Cavidade Óptica.

  Já o segundo, "teste cruzado", testa tal calibração, de modo que a tanto para as calibração dos picos do MZI pares e para os ímpares o      resultado do FSR precisa ser o mesmo.

# Como analisar os seus próprios dados experimentais:

  Este repositório vem com um arquivo de exemplo (dados_mzi.csv), mas o algoritmo foi construído de forma que você possa utilizá-lo no seu    próprio setup! Então, renomeie o seu arquivo para "dados_mzi.csv" para facilitar seu processo.

# Para analisar as suas medições:

  Exporte os dados garantindo que as colunas representem a Transmissão da Cavidade, o sinal do Mach-Zehnder e o espectro do HCN.

  Salve o seu arquivo no formato .csv e renomeie-o exatamente para "dados_mzi.csv".

  Coloque esse novo arquivo dentro da mesma pasta do script apenasgraficos.py ou testecruzado.py, substituindo o arquivo de exemplo.

  Rode o código e acompanhe a geração automática dos gráficos de calibração e do ajuste não-linear da ressonância.

  # Como preparar o ambiente e rodar o código:

  Antes de executar o script de calibração, você precisará instalar as bibliotecas matemáticas e de visualização de dados do Python. Abra o   seu terminal e digite o seguinte comando:

              ```bash
                  pip install numpy pandas matplotlib scipy

  Dependendo da configuração do seu sistema opercional, pode ser necessário utilizar pip3 no lugar de pip.

  E então executar no terminal:

                python apenasgraficos.py 
     ou 
                python testecruzado.py
