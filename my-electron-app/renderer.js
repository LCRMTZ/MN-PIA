const methodOptions = {
    interpolation: [
        { value: "lagrange", text: "Lagrange" },
        { value: "newton-divided", text: "Newton con Diferencias Divididas" },
        { value: "newton-forward", text: "Newton hacia Adelante" },
        { value: "newton-backward", text: "Newton hacia Atrás" }
    ],
    "linear-equations": [
        { value: "montante", text: "Montante" },
        { value: "gauss-seidel", text: "Gauss-Seidel" },
        { value: "gauss-jordan", text: "Gauss-Jordan" },
        { value: "jacobi", text: "Jacobi" },
        { value: "eliminacion-gaussiana", text: "Eliminación Gaussiana" }
    ],
    "nonlinear-equations": [
        { value: "bisection", text: "Bisección" },
        { value: "falsa-posicion", text: "Falsa Posición" },
        { value: "newton-raphson", text: "Newton-Raphson" },
        { value: "punto-fijo", text: "Punto Fijo" },
        { value: "secante", text: "Secante" }
    ],
    "differential-equations": [
        { value: "euler-modificado", text: "Euler Modificado" },
        { value: "segundo-orden", text: "2do Orden" },
        { value: "tercer-orden", text: "3er Orden" },
        { value: "runge-kutta", text: "Runge-Kutta" }
    ],
    integration: [
        { value: "trapezoidal", text: "Regla Trapezoidal (RTMN)" },
        { value: "simpson-1-3", text: "1/3 de Simpson" },
        { value: "simpson-3-8", text: "3/8 de Simpson" },
        { value: "newton-cotes-cerradas", text: "Newton-Cotes Cerradas" }
    ],
    "minimos-cuadrados": [
        { value: "lineal", text: "Mínimos Cuadrados Lineal" },
        { value: "cuadratica", text: "Mínimos Cuadrados Cuadrática" },
        { value: "cubica", text: "Mínimos Cuadrados Cúbica" },
        { value: "linea-recta", text: "Línea Recta" }
    ]
};

const methodTypeSelect = document.getElementById("method-type");
const method1Select = document.getElementById("method1");
const method2Select = document.getElementById("method2");

function populateMethodOptions() {
    const selectedType = methodTypeSelect.value;
    const options = methodOptions[selectedType] || [];

    method1Select.innerHTML = "";
    method2Select.innerHTML = "";

    options.forEach(option => {
        const optionElement1 = document.createElement("option");
        optionElement1.value = option.value;
        optionElement1.text = option.text;
        method1Select.add(optionElement1);

        const optionElement2 = document.createElement("option");
        optionElement2.value = option.value;
        optionElement2.text = option.text;
        method2Select.add(optionElement2);
    });
}

methodTypeSelect.addEventListener("change", () => {
    populateMethodOptions();
    document.getElementById("output").textContent = ""; // Borra resultados previos
    document.getElementById("iteration-progress").innerHTML = ""; // Borra animaciones previas
});

populateMethodOptions(); // Inicializar con el primer tipo de método

// Función para validar los datos de entrada de x y y
function validateInterpolationData(xData, yData) {
    const x = xData.split(',').map(Number);
    const y = yData.split(',').map(Number);

    if (x.includes(NaN) || y.includes(NaN)) {
        return false;
    }

    if (x.length !== y.length) {
        return false;
    }

    return true;
}

// Comparación de métodos
document.getElementById('compare-button').addEventListener('click', async () => {
    try {
        const method1 = document.getElementById('method1').value;
        const method2 = document.getElementById('method2').value;
        const xValue = parseFloat(document.getElementById('x-value').value);
        const xData = document.getElementById('x-data').value;
        const yData = document.getElementById('y-data').value;

        if (isNaN(xValue)) {
            document.getElementById('output').textContent = "Por favor, ingresa un valor válido para x.";
            return;
        }

        if (!validateInterpolationData(xData, yData)) {
            document.getElementById('output').textContent = "Por favor, ingresa puntos válidos para x y y.";
            return;
        }

        const xValues = xData.split(',').map(Number);
        const yValues = yData.split(',').map(Number);

        if (method1 === method2) {
            document.getElementById('output').textContent = "Por favor, elige dos métodos diferentes para comparar.";
            return;
        }

        const [module1, module2] = await Promise.all([
            import(`./methods/${method1}.js`),
            import(`./methods/${method2}.js`)
        ]);

        const result1 = module1.calculateMethod(xValue, xValues, yValues);
        const result2 = module2.calculateMethod(xValue, xValues, yValues);

        document.getElementById('output').innerHTML = `
            <p><strong>${method1}:</strong> Resultado: ${result1.value}, Iteraciones: ${result1.iterations}</p>
            <p><strong>${method2}:</strong> Resultado: ${result2.value}, Iteraciones: ${result2.iterations}</p>
            <p><strong>Comparación:</strong> ${result1.iterations < result2.iterations ? method1 : method2} fue más rápido en cuanto a iteraciones.</p>
        `;

        animateIterations(result1.iterations, result2.iterations);
    } catch (error) {
        console.error(error);
        document.getElementById('output').textContent = "Hubo un error al realizar la comparación.";
    }
});

// Función para animar el progreso de las iteraciones
function animateIterations(iter1, iter2) {
    const method1Progress = document.getElementById('method1-progress');
    const method2Progress = document.getElementById('method2-progress');
    
    method1Progress.style.width = `${(iter1 * 100) / Math.max(iter1, iter2)}%`;
    method2Progress.style.width = `${(iter2 * 100) / Math.max(iter1, iter2)}%`;
}
