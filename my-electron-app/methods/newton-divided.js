export function calculateMethod(xValue, xData, yData) {
    const n = xData.length;
    let result = yData[0];
    let coef = [...yData];

    // Calculamos las diferencias divididas
    for (let j = 1; j < n; j++) {
        for (let i = n - 1; i >= j; i--) {
            coef[i] = (coef[i] - coef[i - 1]) / (xData[i] - xData[i - j]);
        }
    }

    // Evaluamos el polinomio de Newton
    let productTerm = 1;
    for (let j = 1; j < n; j++) {
        productTerm *= (xValue - xData[j - 1]);
        result += coef[j] * productTerm;
    }

    const iterations = (n * (n - 1)) / 2; // Aproximación de número de iteraciones
    return { value: result, iterations: iterations };
}