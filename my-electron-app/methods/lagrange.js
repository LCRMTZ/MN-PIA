// lagrange.js
export function calculateMethod(xValue, xData, yData) {
    let iterations = 0;
    let result = 0;
    const n = xData.length;

    // Cálculo de la interpolación de Lagrange
    for (let i = 0; i < n; i++) {
        let term = yData[i];
        for (let j = 0; j < n; j++) {
            if (j !== i) {
                term *= (xValue - xData[j]) / (xData[i] - xData[j]);
            }
        }
        result += term;
        iterations++;
    }

    return { value: result, iterations };
}
