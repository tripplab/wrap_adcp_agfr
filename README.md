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





