// newton-divided.js
export function calculateMethod(xValue, xData, yData) {
    let iterations = 0;
    const n = xData.length;
    let result = yData[0];
    let diff = [...yData]; // Copia de los valores yData para las diferencias divididas

    // Calcular las diferencias divididas
    for (let i = 1; i < n; i++) {
        for (let j = n - 1; j >= i; j--) {
            diff[j] = (diff[j] - diff[j - 1]) / (xData[j] - xData[j - i]);
        }
    }

    // Usar las diferencias divididas para interpolar
    for (let i = 1; i < n; i++) {
        let term = diff[i];
        for (let j = 0; j < i; j++) {
            term *= (xValue - xData[j]);
        }
        result += term;
        iterations++;
    }

    return { value: result, iterations };
}
