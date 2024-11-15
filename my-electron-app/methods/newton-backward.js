// newton-backward.js
export function calculateMethod(xValue, xData, yData) {
    let iterations = 0;
    let result = yData[yData.length - 1];
    const n = xData.length;

    // Calcular las diferencias divididas hacia atrás
    let diff = [...yData]; // Copia de los valores yData
    for (let i = 1; i < n; i++) {
        for (let j = n - 1; j >= i; j--) {
            diff[j] = diff[j] - diff[j - 1];
        }
    }

    // Interpolar usando las diferencias hacia atrás
    let h = xData[1] - xData[0];  // Asumimos que los intervalos son uniformes
    let p = (xValue - xData[n - 1]) / h;

    for (let i = 1; i < n; i++) {
        let term = diff[n - 1];
        for (let j = 1; j <= i; j++) {
            term *= (p + j - 1);
        }
        result += term;
        iterations++;
    }

    return { value: result, iterations };
}
