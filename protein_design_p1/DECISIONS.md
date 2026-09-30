# Decisiones del dataset — P1: nanobodies unidos a proteína (HER2 marcado)

Bitácora de curación. Cada vez que cambies un umbral o el alcance, apunta aquí
qué cambiaste y por qué, con fecha (apartado 5).

## 1. Alcance: qué entra en el dataset

- **Objeto de estudio**: TODOS los nanobodies (dominios VHH) del PDB que están unidos a
  una proteína, tomados de SAbDab (base de datos de Oxford con las estructuras de
  anticuerpos del PDB). Los que se unen a HER2/ERBB2 se **marcan** con una columna
  (`is_her2 = True`), pero no se separan del resto.
- **Por qué todos y no solo HER2**: en el PDB solo hay 3 nanobodies anti-HER2 con
  calidad utilizable. Con 3 secuencias no se puede deduplicar ni describir nada con
  sentido. Con todos los nanobodies (~900 secuencias distintas) sí, y HER2 sigue
  identificado para reutilizarlo en P3 y P4.
- **Qué es un nanobody, en esta tabla**: SAbDab lo etiqueta como `type = SD-H`
  (single domain, heavy): un anticuerpo con cadena pesada y **sin cadena ligera**. Un
  anticuerpo convencional, como trastuzumab (Herceptin), tiene las dos cadenas.
  `SD-L` y `VNAR` (de tiburón) no son nanobodies de camélido y no entran.
- **Reutilización futura**: base de P2 (clustering para partir en train/test), P3
  (características de la interfaz nanobody-antígeno) y P4 (comparar binders diseñados
  con los reales). El subconjunto HER2 queda siempre localizable por la columna.
- **Fuente de la tabla**: descarga manual de "all structures" de SAbDab
  (`data/raw/sabdab_nano_all.tsv`). La descarga automática no funciona desde que la web
  se rediseñó (devuelve una página web en vez de la tabla).

## 2. Filtros (`01` y `03`) — en el orden en que se aplican

| Paso | Filtro | Umbral | Por qué |
|---|---|---|---|
| 01 | Es nanobody | `type = SD-H` (una sola cadena pesada) | La tabla mezcla todos los anticuerpos; sin esto entran Fabs convencionales |
| 01 | En complejo con antígeno proteico | hay cadena de antígeno y `antigen_type = PROTEIN` | Interesa la interfaz nanobody-proteína, no nanobodies sueltos ni unidos a azúcares/iones |
| 01 | Fuera falsos nanobodies | se quita la fila si su CDR-H3 es idéntica a la de un anticuerpo convencional de la misma tabla | Son cadenas pesadas sueltas de anticuerpos normales mal etiquetadas (p. ej. 8ffj = cadena pesada de trastuzumab) |
| 01 | Marcar HER2 | el nombre del antígeno contiene her2 / erbb2 / erbb-2 → `is_her2 = True` | No filtra; solo etiqueta para P3/P4 |
| 03 | Resolución | ≤ 3.0 Å | Más laxo que el 2.5 Å habitual porque parte de las estructuras nuevas son cryo-EM, que da resoluciones algo peores pero fiables en interfaces |
| 03 | Método experimental | rayos X o cryo-EM | NMR no da una resolución comparable |
| 03 | Longitud del dominio variable (`VH`) | 100–140 aa | Un VHH típico mide unos 110–130 aa; fuera de ese rango suele ser un fragmento o algo mal etiquetado |

Las estructuras de las que RCSB no da resolución o método (por ejemplo NMR) se
descartan y el script 03 avisa de cuántas son. El `filter_log.csv` lleva una columna
`n_her2` para ver cuántos nanobodies anti-HER2 sobreviven a cada paso.

Para probar otro umbral sin tocar el código: `python 03_quality_filter.py --max-res 2.5`.
El nombre de cada línea de `filter_log.csv` incluye el umbral usado.

## 3. Deduplicación (`04`)

- **Herramienta**: script propio en Python (sin instalar nada). Sustituye a CD-HIT,
  que es difícil de instalar en Windows.
- **Identidad**: se alinean las dos secuencias y se cuentan las posiciones con el mismo
  aminoácido, divididas entre la longitud de la más corta (misma definición que CD-HIT).
- **Velocidad**: primero se agrupan las secuencias VH exactamente iguales (una sola vez
  cada una) y luego, antes de cada alineamiento, se descarta el par si comparten muy
  pocos "trozos" de 4 aminoácidos (dos secuencias con ≥95 % de identidad comparten casi
  todos). Es solo un atajo de velocidad: se comprobó con los datos reales que no descarta
  ningún par que llegase al umbral, así que el resultado es idéntico a comparar todo con
  todo. El paso completo tarda unos segundos.
- **Umbral**: 95 % de identidad.
  - 100 % solo quitaría copias exactas y dejaría variantes de 1–2 mutaciones.
  - Un umbral muy bajo (30–40 %) agruparía nanobodies realmente distintos.
  - 95 % quita lo casi idéntico sin borrar diversidad real.
- **Representante de cada grupo**: la estructura con mejor resolución.
- **Limitación conocida**: el armazón de los nanobodies es muy parecido entre sí y lo que
  decide a qué se pegan es la CDR-H3. Comparar la secuencia completa puede juntar en un
  grupo dos nanobodies con CDR-H3 distinta. El script avisa cuando pasa. Mejora posible
  para P3/P4: agrupar solo por CDR-H3 (ya viene anotada en la tabla).
- Para probar otro umbral: `python 04_dedup_identity.py --identity 0.90`.

## 4. Versionado (`05`)

- Cada versión vive en `data/releases/<versión>/` con las secuencias, la tabla completa
  de representantes, la asignación a grupos, el log de filtros y un `MANIFEST.txt`.
- El `MANIFEST.txt` guarda una huella SHA-256 de los archivos: si dentro de un año la
  huella coincide, el archivo es exactamente el mismo.
- El script se niega a sobrescribir una versión existente; para cambiar algo, crea la
  siguiente (v1.1). `CHANGELOG.md` acumula el historial.
- Carpetas + huella en vez de DVC/git-lfs: el dataset son unos pocos KB.

## 5. Historial de cambios de este documento

- **[sept 2026]** Versión inicial. Objetivo cambiado de Lysozyme C a nanobodies anti-HER2 para
  conectar con diseño de anticuerpos (P3/P4).
- **[29 sept 2026]** *Escribe aquí tu línea sobre el cambio de alcance: qué encontraste
  (los primeros 32 resultados eran anticuerpos convencionales, no nanobodies), por qué
  pasó y qué cambia ahora en los filtros.*
