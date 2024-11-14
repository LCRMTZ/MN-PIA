document.getElementById('generate-fields-button').addEventListener('click', () => {
    const nFxValues = parseInt(document.getElementById('n-fx-values').value);
    const nXValues = parseInt(document.getElementById('n-x-values').value);
    let fxText = "";
    let xText = "";

    if (isNaN(nFxValues) || nFxValues <= 0) {
        document.getElementById('funciones').innerHTML = "<p>Error: Ingresa un número válido para los valores fx.</p>";
        return;
    }

    if (isNaN(nXValues) || nXValues <= 0) {
        document.getElementById('funciones').innerHTML = "<p>Error: Ingresa un número válido para los valores x.</p>";
        return;
    }

    // Generar campos para fx
    fxText += "<h3>Valores de fx:</h3>";
    for (let i = 1; i <= nFxValues; i++) {
        fxText += `<p>fx ${i}:</p><input type="number" id="fx-value-${i}" required><br>`;
    }

    // Generar campos para x
    xText += "<h3>Valores de x:</h3>";
    for (let i = 1; i <= nXValues; i++) {
        xText += `<p>x ${i}:</p><input type="number" id="x-value-${i}" required><br>`;
    }

    document.getElementById("funciones").innerHTML = fxText + xText;
});

document.getElementById('compare-button').addEventListener('click', () => {
    const method1 = document.getElementById('method1').value;
    const method2 = document.getElementById('method2').value;

    // Obtener el valor de x para la comparación
    const xValue = parseFloat(document.getElementById('x-value-1').value);

    if (method1 === method2) {
        document.getElementById('output').textContent = "Por favor, elige dos métodos diferentes para comparar.";
        return;
    }

    Promise.all([
        import(`./methods/${method1}.js`),
        import(`./methods/${method2}.js`)
    ]).then(([module1, module2]) => {
        const { calculateMethod: calculateMethod1 } = module1;
        const { calculateMethod: calculateMethod2 } = module2;

        // Obtener el resultado y el número de iteraciones
        const result1 = calculateMethod1(xValue);
        const result2 = calculateMethod2(xValue);

        // Mostrar resultados numéricos
        document.getElementById('output').innerHTML = `
            <p><strong>${method1}:</strong> Resultado: ${result1.value}, Iteraciones: ${result1.iterations}</p>
            <p><strong>${method2}:</strong> Resultado: ${result2.value}, Iteraciones: ${result2.iterations}</p>
            <p><strong>Comparación:</strong> ${result1.iterations < result2.iterations ? method1 : method2} fue más rápido en cuanto a iteraciones.</p>
        `;

        // Crear animación de progreso de iteraciones
        animateIterations(result1.iterations, result2.iterations);
    }).catch(error => {
        console.error("Error al cargar los métodos:", error);
        document.getElementById('output').textContent = "Error al cargar los métodos. Verifique los archivos de métodos.";
    });
});

function animateIterations(iterations1, iterations2) {
    const method1Progress = document.getElementById('method1-progress');
    const method2Progress = document.getElementById('method2-progress');

    // Limpiar cualquier punto de progreso previo
    method1Progress.innerHTML = '';
    method2Progress.innerHTML = '';

    // Crear puntos de progreso para cada método según la cantidad de iteraciones
    for (let i = 0; i < Math.max(iterations1, iterations2); i++) {
        const dot1 = document.createElement('div');
        dot1.classList.add('dot');
        method1Progress.appendChild(dot1);

        const dot2 = document.createElement('div');
        dot2.classList.add('dot');
        method2Progress.appendChild(dot2);
    }

    // Animar el progreso de las iteraciones
    let currentIteration = 0;
    const interval = setInterval(() => {
        if (currentIteration < iterations1) {
            method1Progress.children[currentIteration].classList.add('active');
        }
        if (currentIteration < iterations2) {
            method2Progress.children[currentIteration].classList.add('active');
        }

        currentIteration++;

        // Finalizar la animación cuando se alcanza el máximo de iteraciones de ambos métodos
        if (currentIteration >= Math.max(iterations1, iterations2)) {
            clearInterval(interval);
        }
    }, 500); // Ajusta el intervalo para la velocidad de animación
}
