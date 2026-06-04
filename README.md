# Simulación numérica de abastecimiento, precios y conflicto social en contexto de crisis

Proyecto final de Métodos Numéricos desarrollado como una página web interactiva con HTML, CSS, JavaScript y Chart.js.

## Escenarios incluidos

El proyecto incluye los 7 escenarios solicitados:

1. **Escenario A:** Optimización del abastecimiento y red de transporte mediante sistemas de ecuaciones lineales.
2. **Escenario B:** Vaciado crítico de reservas de carburantes mediante ecuaciones diferenciales ordinarias.
3. **Escenario C:** Curva continua de precios de alimentos mediante interpolación.
4. **Escenario D:** Costo acumulado y pérdida de poder adquisitivo mediante integración numérica.
5. **Escenario E:** Umbrales críticos mediante raíces de ecuaciones.
6. **Escenario F:** Rumores de desabastecimiento y sensibilidad de sistemas mal condicionados.
7. **Escenario G:** Difusión de opinión o descontento social mediante ecuaciones diferenciales.

## Métodos implementados

- Sistemas lineales: LU, Jacobi, Gauss-Seidel, SOR y Gradiente Conjugado.
- Raíces: Bisección, Newton-Raphson y Secante.
- Interpolación: Lagrange, Newton y Splines cúbicos naturales.
- Integración: Trapecio, Simpson 1/3 y Simpson 3/8.
- Ecuaciones diferenciales: Euler, Heun y Runge-Kutta de cuarto orden (RK4).
- Análisis de sensibilidad: número de condición y perturbación de datos.

## Estructura del proyecto

```text
proyecto_metodos_numericos/
├── index.html
├── assets/
│   ├── styles.css
│   └── script.js
├── README.md
├── entrega_proyecto.txt
└── autoevaluacion_individual.txt
```

## Cómo ejecutar localmente

1. Descargar o clonar el repositorio.
2. Abrir `index.html` en un navegador moderno.
3. Editar los datos de los formularios y presionar los botones de simulación.

## Cómo publicar en GitHub Pages

1. Crear un repositorio en GitHub, por ejemplo: `proyecto-metodos-numericos`.
2. Subir todos los archivos del proyecto.
3. Ir a **Settings → Pages**.
4. En **Source**, seleccionar `Deploy from a branch`.
5. Elegir rama `main` y carpeta `/root`.
6. Guardar.
7. GitHub entregará un enlace similar a:

```text
https://usuario.github.io/proyecto-metodos-numericos/
```

## Cómo publicar en Netlify

1. Ingresar a Netlify.
2. Seleccionar **Add new site → Deploy manually**.
3. Arrastrar la carpeta del proyecto.
4. Copiar el enlace generado.

## Nota académica

Los datos iniciales son referenciales y sirven para demostrar la aplicación de los métodos numéricos. Para un estudio real, se deben reemplazar por datos levantados en campo o fuentes oficiales.
