# Identificação de Variantes de Mielofibrose

(Todos os arquivos estão disponíveis em https://github.com/Gnoogit/somatico_2024.git.)

## Alunos responsáveis
- Caio Elias
- Clara Gusmão
- Gabriel Oliveira
- Gilson Borges
- José Filho
- Kaira Tomaz

# Workflow

1. Clonagem do repositório do projeto
2. Instalação dos recursos necessários
3. Filtragem das variáveis de interesse
4. Análise de resultados
5. Conclusão

## 1. Clonagem do repositório do projeto
Com o comando git clone, podemos importar todos os arquivos de uma única vez. O repositório contém as amostras necessárias para análise, assim como scripts extras e arquivos de referência.

```
%%bash
git clone https://github.com/Gnoogit/somatico_2024.git
```

## 2. Instalação dos recursos necessários

a) Instalar bcftools com o plugin split-vep.

O plugin permite extrair os campos de anotações estruturadas como INFO/CSQ criadas por bcftools/csq ou VEP (em nosso caso VEP).

Mais informações: https://samtools.github.io/bcftools/howtos/plugin.split-vep.html

```
%%bash
git clone --recurse-submodules https://github.com/samtools/htslib.git
git clone https://github.com/samtools/bcftools.git
cd bcftools
make
make install
```

b) Instalar udocker.

Udocker é uma ferramenta básica para executar containers docker simples em sistemas sem privilégios de root. Esta abordagem de rodar como sudo será apenas utilizada no workflow pois ele foi elaborado para ser utilizado com o Google Colab e não é recomendada.

Em nosso caso, sempre que utilizarmos o comando udocker rodamos com a opção `docker --allow-root`.

Mais informações: https://indigo-dc.github.io/udocker/

```
%%bash
pip install udocker
udocker --allow-root install
```

c) Download imagem do ensembl-vep.

Ensembl vep é um conjunto de ferramentas para predição de impactos de variantes. Neste workflow, usaremos o comando de filtragem do vep, para filtrar as variáveis de interesse. Como o udocker foi instalado, é possível fazer o download da imagem do vep usando `udocker --allow-root pull`.

Mais informações: https://grch37.ensembl.org/info/docs/tools/vep/index.html

```
%%bash
udocker --allow-root pull ensemblorg/ensembl-vep
```

# 3. Filtragem das variáveis de interesse

A filtragem das variáveis é executada no script vep-gc.sh. O arquivo Myelofibrosis.txt possui uma lista dos genes de interesse. Para mais informações sobre o script, favor consultar o script completo no repositório.

(Neste workflow, foram utilizadas as amostras do projeto LMA Brasil. Os arquivos VCF do projeto foram convertidos previamente da versão do genoma hg19 para hg38 utilizando o programa gatk LiftoverVcf com as posições hg19ToHg38.over.chain da UCSC.)

```
%%bash
sh somatico_2024/vep-gc.sh Myelofibrosis.txt
```

# 4. Análise de resultados

a) Utilizamos então a ferramenta pandas no código abaixo para gerar uma tabela unificada com as variáveis de interesse a partir das saídas do script vep-gc.sh. Pandas é uma ferramenta para análise e manipulação de dados em Python.

Mais informações: https://pandas.pydata.org/

```
import pandas as pd
import os

pd.set_option('display.max_columns', None)
pd.set_option('display.max_colwidth', None)

output_folder = '/content/somatico_2024/my_samples/vep_output/'
output_path = '/content/df_merged.csv'

dataframes = []

for final_outputs in os.listdir(output_folder):
    if final_outputs.endswith('.tsv'):

        output_pathway = os.path.join(output_folder, final_outputs)
        df = pd.read_csv(output_pathway, sep='\t', index_col=False)
        df['Sample'] = final_outputs.split('_')[1]
        dataframes.append(df)

df_merged = pd.concat(dataframes, ignore_index=True)
df_merged.to_csv(output_path, index=False)
```

b) Com a tabela unificada, é possível gerar um gráfico que mostra o número de variantes que cada amostra apresenta. As amostras que não aparecem no gráfico não possuem as variantes de interesse.

```
df_merged_sample_chart = df_merged.value_counts("Sample")
df_merged_sample_chart.plot.bar(x='Sample')
```

c) E também, é possível gerar um gráfico que mostra a distribuição das variantes de interesse por gene. Alguns genes não aparecem no gráfico pois não foram encontradas amostras com as variantes de interesse nos mesmos dentre as amostras analisadas.

```
df_merged_chr_chart = df_merged.value_counts("SYMBOL")
df_merged_chr_chart.plot.pie(y='SYMBOL', figsize=(5, 5), autopct='%1.1f%%', startangle=90)
```

# 5. Conclusão

Assim, de acordo com as informações obtidas através do nosso workflow, podemos demonstrar que dentre as amostras analisadas, apenas 22 das 30 amostras possuem variantes nos genes de interesse, com a maior parte das variantes presentes nos genes JAK2 e CALR.
