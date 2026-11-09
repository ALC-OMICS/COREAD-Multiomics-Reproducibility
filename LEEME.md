# Reproducibilidad computacional COREAD

Versión 0.2.1: código, entradas de medición procesadas, referencias numéricas y
validaciones técnicas. No se distribuyen el manuscrito, figuras finales, leyendas,
tablas suplementarias editoriales ni formularios de envío.

Consulta [README.md](README.md) para instalar el entorno, ejecutar los siete
módulos y comprobar los resultados. Los archivos en `source_data/` y las
referencias CSV/NPZ son datos técnicos necesarios para reproducir los cálculos.
Los cinco Excel OncoPrint del módulo 02 son entradas originales de cBioPortal,
no tablas suplementarias del artículo.

```bash
python -m pip install -r requirements.txt
python scripts/verify_repository.py
python scripts/run_tests.py
python run_all.py --output-dir run_outputs/full
```

Las salidas nuevas se guardan fuera de los archivos fuente y no deben subirse a
GitHub. Un resultado PASS corresponde al modo y entorno realmente probados; no
constituye validación biológica independiente. Los registros de 0.2.0 se conservan
separados de la validación de esta versión en [validation](validation).
