export function calculateMethod(x) {
    // Datos de ejemplo y cálculo simple
    const xValues = [1.7, 2.4, 3.1];
    const yValues = [0.35, 0.87, 1.03];

    let iterations = 0; // Contador de iteraciones
    const h = xValues[1] - xValues[0];
    const s = (x - xValues[2]) / h;

    const deltaY1 = yValues[2] - yValues[1];
    const deltaY2 = yValues[1] - yValues[0];
    const delta2Y1 = deltaY1 - deltaY2;

    const result = yValues[2] + s * deltaY1 + (s * (s + 1) * delta2Y1) / 2;
    iterations += 3; // Ejemplo de conteo de iteraciones

    return { value: result, iterations };
}
