# Proyecto Final Análisis Numérico <br> Universidad Nacional de Colombia <br> 2021-I
## Introducción al modelado de Epidemias en el contexto de la Salud Pública 
### Juan Cañas, José Díaz, David Ortega, Paula Rodríguez, Edward Soto y Fabio Velandia


En este repositorio se junta todo el trabajo realizado y se dividió en cuatro partes, cada una se encuentra en las siguientes carpetas:

<ul>
  <li>Articulo</li>Se encuentran los archivos del articulo final del proyecto.
    <ul>
      <li>Proyecto_An_lisis_Num_rico: Archivo pdf del artículo.
      <li>Tex_Proyecto_An_lisis_Num_rico: Compilado del código .tex del artículo.
    </ul>
  <br><li>Codigo</li>Se encuentra todo el código de las implementaciones de los modelos descritos en el artículo y librerías usadas para generar las diferentes gráficas, simulaciones y resultados. Para ejecutar correctamente los notebooks instalar el archivo requirements.txt. Los diferentes archivos y carpetas son:
    <ul>
      <li>Capeta Img: Contiene las gráficas obtenidas en el código de ejemplo para curvestat y la gráfica de los resultados de la comparación entre SIR y SEI3RD. 
      <li>Carpeta test_data: Contiene los datos de simulación de las curvas de infección para Dinamarca, así como el conjunto de curvas.
      <li>Comparación modelos SIR y SEI3RD: Comparación experimental y teórica de los modelos SIR y SEI3RD. Simulaciones bajo diferentes intervenciones con la población de Bogotá.
      <li>estadisticas_curvas_escenarios_interes: Ejemplo del uso de la librería curvestat para generar estadísticas de curvas y probabilidades de escenarios de interés.
      <li>Modelo SEI3RD Bogota: Construcción del modelo SEI3RD para Bogotá y generación de proyecciones bajo diferentes intervenciones.
      <li>Modelo SIR con control: Construcción del modelo SIR y experimentos relacionados con intervenciones y control.
      <li>SEI3RD_: Simulación del modelo SEI3RD para Bogotá. Este código está basado en la implementación de <a href='https://saludata.saludcapital.gov.co/osb/datos_abiertos_osb/enf-transmisibles/Script_modelo.txt' target='_blank'>Juan Diego Mejía Becerra</a>.
      <li>Vacunación: Implementación de la vacunación al modelo SIR y resultados.
    </ul>
  <br><li>Presentaciones</li>Se encuentran las presentaciones realizadas sobre los temas expuestos, específicamente:
    <ul>
      <li>Estadisticas de curvas: Problemática por estadísticas descriptivas de tiempo fijo y una solución por estadísticas basadas en curvas y probabilidades de escenarios de interés.
      <li>Intro modelos epidemiologicos: Introducción histórica a los modelos epidemiológicos a través del modelo SIR. Exposición de algunos modelos utilizados en Colombia durante la crisis sanitaria de COVID-19.
    </ul>
  <br><li>Propuesta</li>Se encuentran los archivos de la propuesta inicial del proyecto.
    <ul>
      <li>Propuesta_Analisis_Numerico: Archivo pdf de la propuesta.
      <li>Tex_Propuesta Analisis Numerico: Compilado del código .tex de la propuesta.
    </ul>
</ul>
