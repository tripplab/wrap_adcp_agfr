# pepctrl

`pepctrl.py` generates peptide control sequences for a **peptide of interest (POI)**.

It supports three modes:
- **`decoy`**: matched decoys with composition/physicochemical constraints.
- **`scramble`**: permutations of the POI (same residues, reordered).
- **`random`**: random AA20 peptides of the same length.

Author: trippm@tripplab.com [Feb 2026]
---

## Dependencies

`pepctrl.py` uses only Python standard library modules:

- `argparse`
- `csv`
- `json`
- `math`
- `random`
- `re`
- `sys`
- `time`
- `dataclasses`
- `typing`

So there are **no third-party Python package dependencies** for normal use.

### Required runtime

- Python 3.8+ recommended

---

## Installation

### Option 1: run directly (recommended)

Clone the repo and run the script with Python:

```bash
git clone https://github.com/tripplab/pepctrl.git
cd pepctrl
python3 pepctrl.py --help
```

### Option 2: make it executable

```bash
chmod +x pepctrl.py
./pepctrl.py --help
```

---

## Run in a CLI terminal

From the project directory:

```bash
# Show all options
python3 pepctrl.py --help

# Generate 10 decoys and print CSV to terminal
python3 pepctrl.py --poi ACDEFGHIKLMNPQRSTVWY --type decoy --n 10

# Write output file instead of terminal
python3 pepctrl.py --poi ACDEFGHIKLMNPQRSTVWY --type random --n 25 --out random.csv
```

---

## Testing

There is no separate test suite in this repository right now. A practical CLI smoke test is:

```bash
# 1) Verify CLI parsing/help
python3 pepctrl.py --help

# 2) Verify generation/output shape
python3 pepctrl.py --poi ACDEFGHIKLMNPQRSTVWY --type scramble --n 3 --seed 1 --format csv | head
```

If those commands succeed and output looks correct, your local setup is working.

---

## Quick start

```bash
# 1) Generate 25 matched decoys (CSV to stdout)
python3 pepctrl.py --poi ACDEFGHIKLMNPQRSTVWY --type decoy --n 25

# 2) Generate 20 scrambles and write TSV
python3 pepctrl.py --poi ACDEFGHIKLMNPQRSTVWY --type scramble --n 20 --format tsv --out scrambles.tsv

# 3) Generate 50 random controls with fixed seed (reproducible)
python3 pepctrl.py --poi ACDEFGHIKLMNPQRSTVWY --type random --n 50 --seed 42 --out random.csv
```

---

## Input rules

- Provide the POI via `--poi` **or** `--poi-file` (first non-empty line is used).
- Allowed residues are the standard AA20 letters: `ACDEFGHIKLMNPQRSTVWY`.
- `--type` and `--n` are required.

---

## Modes

### 1) `--type decoy`

Best when you want **fair controls** that are not simple shuffles.

Decoy generation can match:
- mean hydrophobicity (`--tol-hydro-mean`),
- aromatic content (`--match-arom`, `--tol-arom-frac`),
- special residues C/P/G (`--match-special`),
- charge-shape and polarization (`--charge-shape`, `--tol-charge-shape`, `--tol-charge-polarization`),
- identity/complexity constraints (`--max-identity`, `--max-run`, `--min-shannon`, motif exclusions, etc.).

Generation strategy is configurable with `--method` (`mutate`, `compose`, `hybrid`) plus tuning flags like `--batch-size`, `--max-attempts`, and `--refine-steps`.

### 2) `--type scramble`

Generates sequence permutations of the POI (same composition, same length).

By default, this mode does **not** apply full decoy filtering. Enable optional filters with:

```bash
--filter-sr
```

### 3) `--type random`

Generates random AA20 peptides with POI-matched length.

Like scramble mode, optional filtering can be enabled with `--filter-sr`.

---

## Output formats

`--format` supports:
- `csv` (default)
- `tsv`
- `fasta`
- `jsonl`

Use `--out -` (default) for stdout, or set a path.

### POI row behavior

- CSV/TSV include a **POI row by default** so POI metrics are visible alongside generated controls.
- Disable with `--no-poi-row`.

---

## Common options

- `--seed <int>`: reproducible outputs.
- `--max-identity-between-decoys <float>`: diversity constraint among returned sequences.
- `--exclude-motif <motif>`: reject sequences containing motif (repeatable).
  - Prefix with `re:` to use regex (example: `--exclude-motif 're:N[^P][ST]'`).
- `--exclude-motif-file <path>`: motif list file, one per line.
- `--allow-partial`: write fewer than `--n` sequences if target cannot be met.
- `-v`, `-vv`, `-vvv`: stderr manifest verbosity.
- `--timing`: include stage timings in manifest output.

---

## Exit codes

- `0`: success, requested output count met.
- `2`: partial result produced (or insufficient candidates without `--allow-partial`).
- `3`: configuration/input error.

---

## Example workflows

### Matched decoys with stricter identity and motifs excluded

```bash
python3 pepctrl.py \
  --poi-file poi.txt \
  --type decoy \
  --n 100 \
  --max-identity 0.25 \
  --exclude-motif RP \
  --exclude-motif 're:N[^P][ST]' \
  --relax \
  --allow-partial \
  --seed 2025 \
  --format csv \
  --out decoys.csv \
  -vv --timing
```

### FASTA output for random controls

```bash
python3 pepctrl.py --poi ACDEFGHIKLMNPQRSTVWY --type random --n 30 --format fasta --out random.fasta
```

---

## Notes

- Charge handling uses a simplified model with configurable histidine partial charge (`--his-charge`) and optional terminal charges (`--include-termini`).
- Hydrophobicity scale currently supports Kyte-Doolittle (`--hydro-scale kd`).





# run_adcp_agfr_replicas_campaign

Author: trippm@tripplab.com [Feb 2026]
---


1) Contexto de aplicación de run_adcp_agfr_replicas_campaign.py

Este script está pensado para correr campañas reproducibles de docking peptídico con AGFR + ADCP dentro de una estructura de experimento con raíz (exp_root) y condiciones organizadas por proteína y péptido.

Estructura de directorios (Mode B, recomendado) 
La convención central es:

<exp_root>/
  inputs/
    receptors/...
    peptides/...
  runs/
    <protein_id>/
      <peptide_id>/
        run_manifest.json
        targets/
          receptor.pdbqt
          receptor.trg
        logs/
          preflight.log
          prepare_receptor.log
          agfr.log
          timing.json
          timing.csv
        replicas/
          replica_001/
            docking/work/...
            logs/adcp.log
            results/...
          replica_002/...
        results/
          README.txt

    Una “condición” = un par (protein_id, peptide_id) (por ejemplo P001/POI).

    Cada condición tiene su propio target per-condition (targets/receptor.trg) y su receptor.pdbqt.

    Los réplicas viven debajo de replicas/replica_###/, cada una con su propio docking/work, logs, results.
Esquema del run_manifest.json
El manifest es la fuente de verdad para describir la condición y resolver rutas:

    campaign: identidad y coherencia (protein_id, peptide_id, ids).

    inputs: receptor_pdb_path, peptide_seq y/o peptide_pdb_path (rutas relativas permitidas).

    paths: run_dir, targets_dir, replicas_dir, results_dir, logs_dir (relativas a exp_root).

    targets: target_trg (también resoluble desde exp_root).

    replicas: planned_count, index_start, suffix_width, seed_base.

    config: parámetros “científicos” de AGFR/ADCP/prepare_receptor.

Reglas/estrategia implícita

    Reproducibilidad por diseño: lo que se corre está descrito en manifest; el script registra ejecución (timing, logs, runtime).

    Separación target vs replicas: el target (.trg) es caro → se comparte por condición; las réplicas sólo hacen ADCP.

    Coherencia estricta de ruta: run_manifest.json debe vivir exactamente en runs/<protein_id>/<peptide_id>/run_manifest.json y coincidir con campaign.{protein_id,peptide_id}.

    Resolución de rutas: toda ruta relativa en el manifest se interpreta contra exp_root.

    Overrides controlados: operacional siempre manda el CLI; lo científico sólo puede diferir del manifest con --allow-cli-overrides.






2) En qué punto se usa el script (cuándo, cómo, desde dónde, qué espera del usuario)

Cuándo usarlo

    Cuando ya tienes definida una condición runs/Pxxx/pep/ con su run_manifest.json, y quieres:

        Validar que todo está bien armado (sin correr nada).

        Inicializar el target (generar receptor.pdbqt y receptor.trg) una sola vez.

        Correr réplicas de ADCP (una o muchas), idealmente con seeds controladas.

Desde dónde ejecutarlo

    Se asume que corres desde algún lugar dentro del experimento, típicamente el exp_root (como hiciste: exp001/), porque:

        El --run-manifest suele pasarse como ruta relativa (runs/P001/POI/run_manifest.json),

        y el script infiere exp_root a partir de la ubicación del manifest.

Ejemplos típicos (desde exp_root):

    Validación:

python3 ../run_adcp_agfr_replicas_campaign.py \
  --run-manifest runs/P001/POI/run_manifest.json \
  --validate-only \
  --verbosity 2

    Inicialización (target/pdbqt) y detenerse:

python3 ../run_adcp_agfr_replicas_campaign.py \
  --run-manifest runs/P001/POI/run_manifest.json \
  --init-only \
  --real-time \
  --verbosity 2

    Ejecutar réplicas:

python3 ../run_adcp_agfr_replicas_campaign.py \
  --run-manifest runs/P001/POI/run_manifest.json \
  --real-time \
  --verbosity 3

Qué espera del usuario

    Que el run_manifest.json exista y sea coherente con su ruta.

    Que las rutas de inputs/ y paths/ sean correctas (o relativas a exp_root).

    Que el ambiente tenga en PATH los binarios: prepare_receptor, agfr, adcp.

    Que entienda la política de overrides:

        si pasa flags científicas en CLI que contradicen el manifest → debe usar --allow-cli-overrides (o corregir el manifest / quitar flags).






3) Salidas que se obtienen (stdout + archivos)

A) STDOUT / mensajes esperables

    Bloque de preflight (siempre): modo, rutas resueltas, rutas de target, flags operacionales, y decisión de gating:

        si va a construir targets o si los va a “SKIP” por existir y ser válidos,

        si está en validate-only/init-only/run,

        warnings (p.ej. --overwrite deprecado → tratado como --overwrite-reps),

        mismatches manifest/CLI cuando aplica.

    Bloques por etapa con timestamps (dependiendo de --verbosity):

        prepare_receptor (si corre),

        agfr (si corre),

        adcp:r### por réplica (si corre),

        y resumen final de éxito/falla.

B) Archivos y artefactos en disco

    Targets por condición (runs/Pxxx/pep/targets/)

    receptor.pdbqt (salida de prepare_receptor)

    receptor.trg (salida de agfr)

    receptor.log u otros auxiliares que AGFR escriba

    Logs por condición (runs/Pxxx/pep/logs/)

    preflight.log (mensajes explícitos de preflight)

    prepare_receptor.log

    agfr.log

    timing.json y timing.csv (con timings por etapa)

    Por réplica (runs/Pxxx/pep/replicas/replica_###/)

    logs/adcp.log (stdout/stderr de ADCP para esa réplica)

    docking/work/… (working directory de ADCP; contiene outputs crudos y temporales)

    results/ (copias/colección “curada”):

        <jobname>_summary.dlg

        <jobname>_out.pdb

        README.txt con inputs + comando + rutas

    Resultados agregados (runs/Pxxx/pep/results/)

    README.txt con resumen de campaña y punteros a replicas/results

    Actualización del run_manifest.json

    Se agrega/actualiza un bloque runtime con:

        last_run_at, status,

        rutas a timing.json/csv,

        listado de replicas, stages, resolved_config,

        y registro de overrides (cuando aplicó --allow-cli-overrides).

Si quieres, te escribo esta síntesis como “README de campaña” listo para pegar en runs/P001/POI/README.md (con comandos concretos y checklist).







4) Ubicación explícita de resultados de docking (poses y energías) para análisis posterior

En Mode B, los resultados “útiles para análisis” quedan por réplica y además el script deja un punto de entrada agregado.

A) Por réplica (lo más directo y recomendado)

Cada réplica escribe y conserva (vía copia) los dos artefactos clave:

Ruta:
runs/<protein_id>/<peptide_id>/replicas/replica_###/results/

Archivos esperados:

    Poses: .../<jobname>_out.pdb
    Contiene las poses dockeadas (las conformaciones/poses finales que produce ADCP).

    Energías / resumen: .../<jobname>_summary.dlg
    Contiene el resumen del docking (scores/energías reportadas por ADCP, ranking, etc.).

Ejemplo concreto (tu caso):

    runs/P001/POI/replicas/replica_001/results/exp001_P001_POI_r001_out.pdb

    runs/P001/POI/replicas/replica_001/results/exp001_P001_POI_r001_summary.dlg

Además, para trazabilidad:

    runs/P001/POI/replicas/replica_###/logs/adcp.log (log completo de ADCP)

    runs/P001/POI/replicas/replica_###/results/README.txt (comando exacto, inputs, rutas)

B) En el working folder (crudo, útil si falta algo en results/)

ADCP genera sus salidas originalmente bajo el workdir de cada réplica:

runs/<protein_id>/<peptide_id>/replicas/replica_###/docking/work/

Ahí típicamente aparecen también:

    <jobname>_summary.dlg

    <jobname>_out.pdb

    (y otros temporales / carpetas internas que ADCP crea)

Si por algún motivo el “copiado a results/” no encontrara un archivo, lo primero es revisar este directorio.

C) Punto de entrada agregado a nivel campaña

El script genera:

runs/<protein_id>/<peptide_id>/results/README.txt

Ahí se listan las rutas de salida de cada réplica (punteros a replica_###/results/), para que puedas ir directo al set completo sin buscar a mano.

D) Para análisis automatizado (globs útiles)

Desde runs/P001/POI/:

    Todas las poses:

find replicas -path "*/results/*_out.pdb" -print

    Todos los resúmenes/energías:

find replicas -path "*/results/*_summary.dlg" -print








# anal_adcp_agfr_replicas_campaign

Author: trippm@tripplab.com [Feb 2026]
---

Análisis reproducible de réplicas ADCP + OpenMM para docking de péptidos
Este paquete/script (anal_adcp_agfr_replicas_campaign) analiza resultados de docking de péptidos producidos con AutoDock CrankPep (ADCP) y re‐ranking por minimización con OpenMM (con solvente implícito y score de interacción -reint). Su objetivo principal es transformar un conjunto de réplicas dispersas en disco en:
Tablas consolidadas (CSV) con métricas por réplica y por grupo (proteína–péptido)
Métricas de discriminación entre péptidos (p.ej. POI vs controles) con significancia estadística
Exports de poses para inspección estructural (winners y top-k concatenado)
El enfoque es de generación de hipótesis / priorización, no un estimador concluyente de afinidad o especificidad. El score de OpenMM se usa como criterio relativo bajo un protocolo constante para POI y controles, y se complementa con señales de convergencia estructural (clusterización por contactos).

1. Filosofía y supuestos
1.1 Qué se intenta medir
El pipeline busca evidencia del tipo:
“Para esta proteína, el POI converge repetidamente a un sitio/modo de contacto dominante y los controles no.”
“Para esta proteína, el POI tiende a puntuar mejor (más negativo) que el fondo de controles bajo el mismo protocolo.”
“En esta proteína, todos los péptidos se comportan igual → no hay señal discriminante.”
1.2 Qué NO se debe concluir
El score de OpenMM tras minimización NO es energía libre de unión (ΔG).
Un score más negativo no prueba afinidad específica por sí solo.
La minimización puede introducir artefactos si se interpreta como termodinámica real.
Por eso, anal_adcp_agfr_replicas_campaign reporta:
distribuciones (medianas/IQR, no “best run”),
tamaños de efecto y pruebas estadísticas,
convergencia y “modo dominante”.

2. Estructura esperada del experimento en disco
anal_adcp_agfr_replicas_campaign está diseñado para recorrer un experimento con estructura tipo:
exp001/
 analysis/                    (salidas de anal_adcp_agfr_replicas_campaign; puede existir o crearse)
 runs/
   P001/
     POI/
       replicas/replica_001/results/*_summary.dlg
       replicas/replica_001/results/*_omm_rescored_out.pdb   (si existe)
       replicas/replica_001/docking/work/*_out.pdb           (opcional)
       ...
     polya/
     random001/
     scramble001/
     decoy001/
     x12/
   P002/
   ...
2.1 Auto-discovery
El script no requiere que existan resultados en todas las carpetas. Recorre todo lo que haya bajo runs/ y:
procesa lo que encuentre (réplicas con *_summary.dlg)
emite warnings resumidos por grupo vacío o incompleto
continúa sin interrumpirse
2.2 Archivos clave por réplica (fuentes de verdad)
Para una carpeta .../replica_XXX/, los archivos relevantes son:
Canónico: results/*_summary.dlg
Contiene:
tabla de clusters ADCP
energías OpenMM por modelo minimizado
tabla final OMM Ranking (re-ranking por interacción)
Opcional: results/*_omm_rescored_out.pdb
PDB multi-modelo reordenado por ranking OpenMM (modelo 1 = rank 1), con USER: incluyendo dE_Interaction.
Opcional: docking/work/*_out.pdb
PDB con modelos ADCP pre-minimización, USER: SCORE ... por modelo.
Si un archivo no existe, anal_adcp_agfr_replicas_campaign utiliza lo disponible y marca QA.

3. Cómo se define el “winner” por réplica
Cada réplica produce varios modelos (poses). anal_adcp_agfr_replicas_campaign define:
Winner post-OpenMM = la fila con RankOpenMM = 1 en la tabla OMM Ranking del *_summary.dlg.
Este winner corresponde, cuando existe *_omm_rescored_out.pdb, a MODEL 1 (porque el PDB está reordenado por el ranking).
3.1 Score principal por réplica
El score primario para análisis es:
winner_omm_dE_interaction
que corresponde a:
dE_Interaction reportado por OpenMM, equivalente a
E_complex − E_receptor − E_peptide
bajo solvente implícito (según el protocolo).
Interpretación:
más negativo → interacción más favorable según este criterio de minimización
se usa comparativamente (POI vs controles) y con réplicas

4. Identificación de clusters y convergencia (definición B)
ADCP reporta una tabla de clusters “mode | affinity | clust size | best run …”.
El OMM Ranking reporta best run asociado al modelo. Para asignar un identificador de cluster al winner:
Definición B (adoptada):
winner_best_run_pose_id se toma del renglón winner en OMM Ranking.
Se busca en la tabla de clusters ADCP el mode cuyo best run coincide.
Ese mode se guarda como winner_cluster_mode_id.
Esto permite calcular convergencia por grupo (proteína–péptido) usando la distribución de winner_cluster_mode_id en las réplicas.

5. Archivos de salida (en analysis/)
anal_adcp_agfr_replicas_campaign produce varias tablas CSV y archivos auxiliares.
5.1 replicas_parsed.csv (tabla principal por réplica)
Una fila = 1 réplica (winner OpenMM)
Campos típicos (pueden variar ligeramente si faltan archivos):
IDs: experiment_id, protein_id, peptide_id, replica_id
Trazabilidad: summary_path, rescored_pdb_path, out_pdb_path
Parámetros: sequence, N_runs, n_evals, nmin, env, nitr
ADCP: adcp_bestEnergy, adcp_bestEnergy_run
Winner:
winner_omm_dE_interaction (score principal)
winner_omm_dE_complex_minus_receptor (auxiliar)
winner_adcp_affinity
winner_rank_adcp
winner_cluster_size
winner_best_run_pose_id
winner_cluster_mode_id (por definición B)
QA: qa_status, qa_message
Cómo usarla:
Para comparar POI vs control: filtrar por protein_id y peptide_id, graficar la distribución de winner_omm_dE_interaction.
Para convergencia: contar frecuencia de winner_cluster_mode_id.

5.2 topk_parsed.csv (top-k re-rankeado por réplica)
Una fila = 1 modelo (pose) dentro del top-k (nmin)
Incluye para cada réplica y rank_openmm=1..k:
omm_dE_interaction, adcp_affinity, cluster_size, best_run_pose_id, cluster_mode_id, etc.
Cómo usarla:
Ver si el POI mantiene ventaja en top-3, no solo en top-1.
QA de reordenamiento: ver diferencias entre rank_openmm y rank_adcp.

5.3 clusters_parsed.csv (tabla de clusters ADCP por réplica)
Una fila = 1 cluster mode (1..10) por réplica
mode, adcp_affinity, cluster_size, best_run_pose_id
Cómo usarla:
Diagnóstico de landscape: ¿hay un modo dominante grande o muchos modos pequeños?
Comparar “estructura del clustering” entre péptidos.

5.4 group_summary.csv (resumen por proteína–péptido)
Una fila = 1 grupo (protein_id, peptide_id)
Se calcula usando replicas_parsed.csv con qa_status == OK.
Incluye:
Tamaños
n_replicas_total, n_ok, fail_rate
Energía (score principal)
Eint_median, Eint_IQR (y opcionales mean/sd/MAD)
Convergencia
cluster_mode_top
cluster_mode_top_freq (= f1_cluster)
cluster_mode_entropy
cluster_mode_Neff (= exp(entropy))
n_unique_modes
Diagnósticos auxiliares
winner_cluster_size_median, winner_cluster_size_IQR
Cómo usarla:
Comparar rápidamente POI vs polya en una proteína:
si POI tiene Eint_median más favorable y
cluster_mode_top_freq más alto / Neff más bajo → señal de convergencia.

5.5 discrimination.csv (comparación entre péptidos dentro de proteína)
Una fila = 1 comparación (A vs B) dentro de una proteína
Usa los winner_omm_dE_interaction de cada grupo:
wins_A_over_B (proporción de pares con A “mejor”)
AUC (equivalente a Mann–Whitney U normalizado)
cliffs_delta (= 2*AUC − 1)
p_mannwhitney
p_permutation (si se activa)
p_adj_fdr (cuando hay múltiples comparaciones)
convergencia: f1_A, f1_B, Neff_A, Neff_B, delta_*
verdict (diagnóstico automático)
Cómo interpretar:
AUC ~ 0.5 → sin discriminación
AUC → 1.0 → A casi siempre mejor que B
cliffs_delta magnitud del efecto (regla práctica):
|δ| < 0.147: pequeño
0.147–0.33: pequeño–medio
0.33–0.474: medio
0.474: grande

6. Exports de poses (PDB)
6.1 Winners por réplica
Directorio:
analysis/poses_winners/<protein_id>/<peptide_id>/
Se exporta 1 PDB por réplica (winner), típicamente:
desde *_omm_rescored_out.pdb → MODEL 1
si falta, anal_adcp_agfr_replicas_campaign puede omitir el export y reportarlo en QA.
Uso:
inspección visual rápida de convergencia de sitio y geometría.
6.2 Top-k concatenado (todas las réplicas)
Directorio:
analysis/poses_topk_concat/<protein_id>/<peptide_id>/
Archivos:
topk_concat.pdb (por grupo proteína–péptido)
Contiene:
si hay 30 réplicas y k=5 → hasta 150 modelos.
Orden recomendado:
réplica_001 rank1..rank5
réplica_002 rank1..rank5
...
Cada MODEL incluye encabezado USER: con energías, para rastreabilidad.
Uso:
generar mapas de contactos, clustering externo, o scripts de visualización que recorran modelos.

7. Filtro de sanity para energías
anal_adcp_agfr_replicas_campaign puede marcar (no necesariamente descartar) réplicas con:
winner_omm_dE_interaction NaN / no parseable
valores fuera de un rango razonable (umbral configurable)
Objetivo:
evitar que 1–2 outliers numéricos arruinen pruebas o mediana/IQR.
mantener trazabilidad (réplica se marca como QA fail con mensaje).

8. Ejecución
8.1 CLI
Ejemplo:
anal_adcp_agfr_replicas_campaign --exp-root exp001 --outdir analysis
Donde:
--exp-root apunta al directorio del experimento (contiene runs/)
--outdir es donde escribir outputs (por convenio analysis/)
8.2 Comportamiento con datos incompletos
Si existen carpetas random001/ etc. sin réplicas parseables:
anal_adcp_agfr_replicas_campaign no falla
reporta algo como:
runs/P001/random001: 0 replicas parsed (no summary.dlg found)
y continúa.

9. Flujo recomendado de uso e interpretación
Correr anal_adcp_agfr_replicas_campaign tras obtener réplicas para POI y un control (p.ej. polya).
Ver analysis/parse_report.* para asegurar:
n_ok ~ número esperado de réplicas
fail_rate bajo
Revisar group_summary.csv:
comparar Eint_median y Neff/f1 entre POI y control
Revisar discrimination.csv:
AUC, Cliff’s delta y p-values
Si hay señal:
escalar a más controles (scramble/random/decoys matched)
y a más proteínas
Si no hay señal:
revisar QA
revisar distribución top-k
inspeccionar poses winners/topk_concat para detectar “pegajosidad” o sitios triviales

10. Notas prácticas y advertencias
Comparaciones entre proteínas distintas deben hacerse con cuidado: las energías de interacción de OpenMM pueden tener offsets dependientes del sistema.
Lo más defendible es comparar dentro de la misma proteína: POI vs controles.
La convergencia por cluster_mode_id depende del clustering interno de ADCP (cutoff, etc.).
Aun así, el contraste POI vs controles bajo el mismo protocolo es informativo.
El pipeline es útil como priorización: para decidir qué proteínas/peptidos merecen rescoring más caro (MM/GBSA, MD corta, etc.) o validación experimental.

11. Archivos auxiliares (recomendados)
Además de los CSV, anal_adcp_agfr_replicas_campaign genera:
analysis/parse_report.txt o .json
Resumen de grupos encontrados, réplicas OK/fallidas y causas.
analysis/run_info.json
Versiones, fecha, argumentos CLI, y parámetros de análisis (k, umbrales QA, etc.).

12. Preguntas frecuentes
¿Qué métrica debo usar para “ranking” final?
Para ranking dentro de una proteína, use una combinación de:
Eint_median (más negativo mejor)
cluster_mode_top_freq alto / Neff bajo (convergencia)
estabilidad (IQR no enorme, fail_rate bajo)
¿Qué pasa si polya “gana” en energía?
Eso puede indicar un régimen de pegajosidad hidrofóbica. No se oculta: se reporta. En ese caso:
los decoys matched se vuelven aún más importantes,
y conviene revisar si el box/entorno favorece superficies hidrofóbicas no específicas.
¿Por qué AUC/Mann–Whitney y no t-test?
Porque:
no asumimos normalidad
el criterio es ordinal (“¿A suele ser mejor que B?”)
es más robusto para distribuciones sesgadas o con colas.




