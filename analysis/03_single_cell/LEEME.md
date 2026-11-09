# Guía rápida

Este repositorio corresponde al análisis **single-cell transcriptomics** de los 466 candidatos funcionales de COREAD. Incluye la data numérica de **Supplementary Table 10** y **Figure 3A–D**, junto con archivos fuente, reglas de priorización, verificaciones y documentación.

El análisis principal usa **pseudobulk a nivel de paciente**. No se trataron las células individuales como réplicas biológicas independientes. La comparación principal es epitelio derivado de tumor frente a epitelio adyacente normal en 36 pacientes emparejados. Los análisis de >=100 células y de composición epitelial son sensibilidades separadas. El análisis MMRd versus MMRp es secundario.

La priorización final conserva los 466 genes y los clasifica en 65 Tier A, 134 Tier B y 267 Tier C. De los 65 Tier A, 47 mantienen enriquecimiento tumoral positivo y FDR <0.05 después de balancear la composición de subtipos epiteliales; estos 47 constituyen el conjunto recomendado para el siguiente análisis espacial.

Para verificar el paquete:

```bash
python -m pip install -r requirements.txt
python run_all.py --verify-only
```

Para regenerar Figure 3A–D desde sus archivos CSV fuente:

```bash
python run_all.py --figures-only --output build-figures
```

El H5 de GSE178341 no se incluye porque pesa ~1.2 GB. El SHA-256 esperado y las instrucciones para reconstruirlo desde partes están en `data/raw/README.md`.
