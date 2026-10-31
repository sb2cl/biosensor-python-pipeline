# Publicar el pipeline Python en GitHub y ejecutarlo desde el navegador

## 1. Repositorio de código: público, pero solo con datos sintéticos

Nombre sugerido: `sb2cl/biosensor-python-pipeline`.

En **GitHub Desktop**:

1. Descomprime `BioSensor_Python_Pipeline_Starter_v0_1.zip`.
2. `File > Add local repository...` y selecciona la carpeta descomprimida.
3. Si informa que no es un repositorio Git, utiliza **Create a repository here**.
4. Realiza un commit inicial: `Initial reproducible Python biosensor pipeline`.
5. Pulsa **Publish repository**, selecciona organización `sb2cl` (si tienes permiso) y decide `Public` para habilitar la demo en MyBinder.
6. Verifica la pestaña **Actions**: debe pasar `Python data pipeline checks`.

**No subas** las hojas Excel originales ni los `.mat` de experimentos o posteriores al repositorio público. La carpeta `private_data/` está excluida por `.gitignore`, pero revisa `git status` antes de publicar.

## 2. GitHub Codespaces: la alternativa oficial de GitHub

1. Abre tu repositorio en GitHub.
2. Haz clic en **Code > Codespaces > Create codespace on main**.
3. Espera a que `postCreateCommand` instale las dependencias.
4. Abre `notebooks/01_excel_to_feature_studio.ipynb` y pulsa **Run All** (kernel Python recomendado).
5. Los CSV aparecen en `outputs/notebook_demo/`.
6. Para usar un Excel real, cárgalo solo en tu sesión privada de Codespaces (o trabaja en un repositorio privado), crea una configuración adaptada y no lo confirmes en Git.

Codespaces tiene recursos y facturación sujetos a la configuración de tu cuenta/organización. Para simulaciones largas puede ser preferible tu ordenador personal o HPC del laboratorio.

## 3. MyBinder: demo pública sin instalar nada

Con el repositorio publicado y `requirements.txt`, abre:

`https://mybinder.org/v2/gh/sb2cl/biosensor-python-pipeline/HEAD?urlpath=lab/tree/notebooks/01_excel_to_feature_studio.ipynb`

MyBinder es gratuito pero **efímero**, con limitaciones de memoria y tiempo, y no es la plataforma apropiada para una inferencia PyVBMC prolongada. Su repositorio debe estar accesible a Binder.

## 4. Primer Excel experimental REAL

Actualmente no hay `.xlsx` originales en `scripts_final.zip`. Para automatizar la lectura real, necesitaré al menos uno de estos archivos o su estructura de hojas, pozos y calibración. Mientras tanto, el ejemplo de la carpeta `examples/` es exclusivamente sintético y el código hace cumplir un mapa de columnas/dosis explícito.

**FdeR:** ventana de observación en **índices MATLAB [72:108]**, inclusivos, sobre series ya procesadas y alineadas. **pLux:** ventana `time_OD` 360-480 minutos. No intercambiarlas.

## 5. Figuras congeladas y modelo

`notebooks/02_frozen_matlab_ppc_figures.ipynb` genera figura con PPC del 21/09/2026 a partir de `.mat` y CSV archivados, sin MATLAB. Para diseño F3/F5 podemos añadir después scripts específicos sobre las matrices y CSV congelados.

La migración del **modelo final de cinco estados, 11 parámetros, dos contextos experimentales** se hará como un módulo separado, con verificación en parámetros idénticos y logp histórico. Solo entonces consideraremos PyBADS/PyVBMC para nuevos ajustes.

## 6. Nuevos modelos ODE en v0.2 (31/10/2026)

Se añaden tres notebooks: `03_minimal_FdeR_RBS.ipynb`,
`04_final_FdeR_MATLAB_port.ipynb` y `05_historical_2025_candidate.ipynb`.
El primero funciona con datos sintéticos y muestra la dirección del cambio de
EC50 al cambiar la RBS del TF. El segundo reproduce el modelo final de 11
parámetros con contexto sintético público, y permite auditoría numérica contra
los `.mat` originales solo si se copian localmente a `private_data/`. El
tercero transcribe el modelo de desrepresión de 2025, con limitaciones
archivísticas explícitas.

**Importante:** La versión de dos RBS del modelo mínimo es una demostración,
no una refit ni una recuperación de las predicciones históricas. Las
observaciones de DR/EC50 reales NO se han añadido al repositorio público.
Para que el modelo simple ajuste el biosensor original necesitaremos los datos
reales y establecer una calibración separada y trazable.

Una vez publicado, MyBinder puede abrir directamente:

`https://mybinder.org/v2/gh/sb2cl/biosensor-python-pipeline/HEAD?urlpath=lab/tree/notebooks/03_minimal_FdeR_RBS.ipynb`

Para experimentar con los `.mat` privados o ejecutar muchas curvas, usar un
**Codespace privado**, ordenador local o HPC, no un Binder público.
