# Workflow - análise de variantes somáticas
---
## Enunciado

Objetivo geral: análise de variantes somáticas detectadas em genes de alto risco
que são fatores de prognósticos adversos em Mielofibrose (MF).

**Identificar as amostras com alterações no TP53 e nos genes de alto risco
( EZH2, CBL, U2AF1, SRSF2, IDH1, IDH2, NRAS ou KRAS ) que são
fatores de prognósticos adversos em mielofibrose.**
---

1. Clonar repositório `clarasgusmao/somatico_2024`
2. Instalar `bcftools +split-vep`
3. Instalar `udocker`
4. Filtrar o VCF com `filter_vep`:

  ```
  -filter "(MAX_AF <= 0.01 or not MAX_AF) and
  (FILTER = PASS or not FILTER matches strand_bias,weak_evidence) and
  (SOMATIC matches 1 or (not SOMATIC and CLIN_SIG matches pathogenic)) and
  (not CLIN_SIG matches benign) and \
  (not IMPACT matches LOW) and \
  (Symbol in hpo/$HPO)
  ```

5. Filtrar Cobertura Total e Frequência Alélica da variante com: `bcftools +split-vep`:
  - `DP>=20 AND AF>=0.1`
6. Resultado: `*.vep.filter.tsv`

The Human Phenotype Ontology (HPO)

1. Myelofibrosis: https://hpo.jax.org/app/browse/term/HP:0011974
2. Abnormal mast cell morphology: https://hpo.jax.org/app/browse/term/HP:0100494
3. Genes de alto risco que são fatores prognósticos adversos em mielofibrose (enunciado)
---
