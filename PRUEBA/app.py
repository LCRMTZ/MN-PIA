from flask import Flask, request
import numpy as np
from scipy.interpolate import lagrange

app = Flask(__name__)

# Verificar si los intervalos son uniformes
def intervalos_uniformes(valores):
    if len(valores) < 2:
        return True
    diferencias = [valores[i+1] - valores[i] for i in range(len(valores) - 1)]
    return all(abs(diferencias[0] - d) < 1e-9 for d in diferencias)


# Interpolación Lineal
def interpolacion_lineal(x, y, valor):
    for i in range(len(x) - 1):
        if x[i] <= valor <= x[i+1]:
            return y[i] + (y[i+1] - y[i]) * (valor - x[i]) / (x[i+1] - x[i])
    return None

# Interpolación de Lagrange
def interpolacion_lagrange(x, y, valor):
    polinomio = lagrange(x, y)
    return polinomio(valor)

# Newton Adelante
def newton_adelante(x, y, valor):
    n = len(x)
    diferencias = [y[:]]
    for i in range(1, n):
        diferencias.append([
            (diferencias[i-1][j+1] - diferencias[i-1][j]) / (x[j+i] - x[j])
            for j in range(n-i)
        ])
    resultado = diferencias[0][0]
    prod = 1
    for i in range(1, n):
        prod *= (valor - x[i-1])
        resultado += prod * diferencias[i][0]
    return resultado

# Newton Atrás
def newton_atras(x, y, valor):
    n = len(x)
    diferencias = [y[:]]
    for i in range(1, n):
        diferencias.append([
            (diferencias[i-1][j] - diferencias[i-1][j-1]) / (x[j] - x[j-i])
            for j in range(i, n)
        ])
    resultado = diferencias[0][-1]
    prod = 1
    for i in range(1, n):
        prod *= (valor - x[-i])
        resultado += prod * diferencias[i][-1]
    return resultado

def newton_diferencias_divididas(x, y, valor):
    n = len(x)
    # Verificar que haya suficientes puntos
    if n < 2:
        raise ValueError("Se necesitan al menos dos puntos para interpolar.")
    
    # Inicializar tabla de diferencias divididas
    coeficientes = [y[:]]  # Copia inicial de los valores de y
    for i in range(1, n):
        fila_actual = []
        for j in range(n - i):
            if (x[j+i] - x[j]) == 0:
                raise ValueError("División por cero en las diferencias divididas.")
            diferencia = (coeficientes[i-1][j+1] - coeficientes[i-1][j]) / (x[j+i] - x[j])
            fila_actual.append(diferencia)
        coeficientes.append(fila_actual)
    
    # Calcular el valor interpolado
    resultado = coeficientes[0][0]
    producto = 1
    for i in range(1, n):
        producto *= (valor - x[i-1])
        resultado += coeficientes[i][0] * producto
    return resultado

# Ruta principal para mostrar el formulario
@app.route('/')
def index():
    with open("index.html", "r") as f:
        html_content = f.read()
    return html_content


# Ruta para procesar datos
@app.route('/process', methods=['POST'])
def process():
    
    try:
        # Leer los datos del formulario
        x_values = list(map(float, request.form['x_values'].split(',')))
        y_values = list(map(float, request.form['y_values'].split(',')))
        test_value = float(request.form['test_value'])
        method_1 = request.form['method_1']
        method_2 = request.form['method_2']

        if len(x_values) != len(y_values):
            return "Error: Las listas X e Y deben tener el mismo tamaño."

        # Verificar uniformidad de los datos
        x_uniforme = intervalos_uniformes(x_values)

        # Diccionario de métodos según la uniformidad
        metodos = {
            "lineal": interpolacion_lineal,
            "lagrange": interpolacion_lagrange,
            "newton_adelante": newton_adelante if x_uniforme else None,
            "newton_atras": newton_atras if x_uniforme else None,
            "newton_diferencias_divididas": newton_diferencias_divididas
        }

        # Resultados de los métodos seleccionados
        resultados = {}
        for metodo in [method_1, method_2]:
            if metodo in metodos and metodos[metodo]:
                resultados[metodo] = metodos[metodo](x_values, y_values, test_value)
            else:
                resultados[metodo] = "No aplicable (intervalos no compatibles)."

        # Preparar salida
        resultado_html = f"""
        <h1>Resultados</h1>
        <ul>
            <li>Método {method_1}: {resultados[method_1]}</li>
            <li>Método {method_2}: {resultados[method_2]}</li>
        </ul>
        <p>Intervalos uniformes: {'Sí' if x_uniforme else 'No'}</p>
        <a href="/">Volver</a>
        """
        return resultado_html
    except Exception as e:
        return f"Error: {str(e)}"


if __name__ == '__main__':
    app.run(debug=True)


