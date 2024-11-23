import tkinter as tk
from tkinter import messagebox
import numpy as np
from scipy.interpolate import lagrange
import sympy as sp

# Verificar si los intervalos son uniformes
def intervalos_uniformes(valores):
    if len(valores) < 2:
        return True
    diferencias = [valores[i+1] - valores[i] for i in range(len(valores) - 1)]
    return all(abs(diferencias[0] - d) < 1e-9 for d in diferencias)


# Interpolación Lineal
# def interpolacion_lineal(x, y, valor):
#     for i in range(len(x) - 1):
#         if x[i] <= valor <= x[i+1]:
#             return y[i] + (y[i+1] - y[i]) * (valor - x[i]) / (x[i+1] - x[i])
#     return None
def interpolacion_lineal(x, y, valor):
    # Validar que x y y tienen más de un valor
    if len(x) <= 1 or len(y) <= 1:
        raise ValueError("x e y deben contener más de un valor cada uno.")
    
    # Validar que x e y tienen la misma cantidad de valores
    if len(x) != len(y):
        raise ValueError("x e y deben tener la misma cantidad de valores.")
    
    # Validar que x no tiene valores repetidos
    if len(set(x)) != len(x):
        raise ValueError("x no debe contener valores repetidos.")
    
    # Convertir a una lista de coordenadas y ordenar por x
    coordenadas = sorted(zip(x, y), key=lambda punto: punto[0])
    
    # Separar nuevamente las listas x e y ordenadas
    x_ordenado, y_ordenado = zip(*coordenadas)
    
  
    # Realizar interpolación lineal
    for i in range(len(x_ordenado) - 1):
        if x_ordenado[i] <= valor <= x_ordenado[i+1]:
            return y_ordenado[i] + (y_ordenado[i+1] - y_ordenado[i]) * (valor - x_ordenado[i]) / (x_ordenado[i+1] - x_ordenado[i])
              
    # Si no se encuentra el intervalo, retornar None
    return None


# Interpolación de Lagrange
def interpolacion_lagrange(x, y, valor):
    # Validar que x y y tienen más de un valor
    if len(x) <= 1 or len(y) <= 1:
        raise ValueError("x e y deben contener más de un valor cada uno.")
    
    # Validar que x e y tienen la misma cantidad de valores
    if len(x) != len(y):
        raise ValueError("x e y deben tener la misma cantidad de valores.")
    
    # Validar que x no tiene valores repetidos
    if len(set(x)) != len(x):
        raise ValueError("x no debe contener valores repetidos.")
    
    polinomio = lagrange(x, y)
    return polinomio(valor)

#Newton hacia adelante
def newton_adelante(x, y, valor):
    # Validar que x y y tienen más de un valor
    if len(x) <= 1 or len(y) <= 1:
        raise ValueError("x e y deben contener más de un valor cada uno.")
    
    # Validar que x e y tienen la misma cantidad de valores
    if len(x) != len(y):
        raise ValueError("x e y deben tener la misma cantidad de valores.")
    
    # Validar que x no tiene valores repetidos
    if len(set(x)) != len(x):
        raise ValueError("x no debe contener valores repetidos.")
    
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
    # Validar que x y y tienen más de un valor
    if len(x) <= 1 or len(y) <= 1:
        raise ValueError("x e y deben contener más de un valor cada uno.")
    
    # Validar que x e y tienen la misma cantidad de valores
    if len(x) != len(y):
        raise ValueError("x e y deben tener la misma cantidad de valores.")
    
    # Validar que x no tiene valores repetidos
    if len(set(x)) != len(x):
        raise ValueError("x no debe contener valores repetidos.")
    
    n = len(x)
    diferencias = [y[:]]
    
    # Construir tabla de diferencias divididas hacia atrás
    for i in range(1, n):
        diferencias.append([
            (diferencias[i-1][j] - diferencias[i-1][j-1]) / (x[j] - x[j-(i+1)])
            for j in range(i, n)
        ])
    
    # Evaluar el polinomio en el valor deseado
    resultado = diferencias[0][-1]
    prod = 1
    for i in range(1, n):
        prod *= (valor - x[-i])
        resultado += prod * diferencias[i][-1]
    
    return resultado


def newton_diferencias_divididas(x, y, valor):
    # Validar que x y y tienen más de un valor
    if len(x) <= 1 or len(y) <= 1:
        raise ValueError("x e y deben contener más de un valor cada uno.")
    
    # Validar que x e y tienen la misma cantidad de valores
    if len(x) != len(y):
        raise ValueError("x e y deben tener la misma cantidad de valores.")
    
    # Validar que x no tiene valores repetidos
    if len(set(x)) != len(x):
        raise ValueError("x no debe contener valores repetidos.")
    
    n = len(x)
    # Verificar que haya suficientes puntos
    if n < 2:
        raise ValueError("Se necesitan al menos dos puntos para interpolar.")
    
    diferencias = [x[i+1] - x[i] for i in range(len(x) - 1)]
    if not all(diff == diferencias[0] for diff in diferencias):
        raise ValueError("Los valores en x no están espaciados de forma uniforme.")

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

# Funciones de resolución de sistemas de ecuaciones lineales
def gauss_jordan(A, b):
    n = len(A)
    M = np.hstack([A, b.reshape(-1, 1)])  # Matriz aumentada [A | b]
    for i in range(n):
        M[i] = M[i] / M[i, i]  # Dividir fila por el pivote
        for j in range(n):
            if j != i:
                M[j] -= M[i] * M[j, i]
    return M[:, -1]

def gauss_eliminacion(A, b):
    n = len(A)
    augmented_matrix = np.hstack([A, b.reshape(-1, 1)])
    for i in range(n):
        max_row = max(range(i, n), key=lambda r: abs(augmented_matrix[r, i]))
        augmented_matrix[[i, max_row]] = augmented_matrix[[max_row, i]]
        for j in range(i + 1, n):
            factor = augmented_matrix[j, i] / augmented_matrix[i, i]
            augmented_matrix[j, i:] -= factor * augmented_matrix[i, i:]

    x = np.zeros(n)
    for i in range(n - 1, -1, -1):
        x[i] = (augmented_matrix[i, -1] - np.dot(augmented_matrix[i, i + 1:n], x[i + 1:])) / augmented_matrix[i, i]
    return x

def Montante(A, b):
    # Convertimos A y b en matrices flotantes para operaciones
    A = np.array(A, dtype=float)
    b = np.array(b, dtype=float).reshape(-1, 1)
    
    # Construimos la matriz ampliada
    Augmented = np.hstack((A, b))
    n = len(A)
    
    # Inicializamos el pivote
    pivot = 1.0

    # Aplicamos el método de Montante
    for k in range(n):
        # Guardamos el pivote actual
        pivot_k = Augmented[k, k]
        
        # Actualizamos toda la matriz
        for i in range(n):
            for j in range(n + 1):
                if i != k and j != k:
                    Augmented[i, j] = (pivot * Augmented[i, j] - Augmented[i, k] * Augmented[k, j]) / pivot_k
        
        # Actualizamos la fila pivote
        for j in range(n + 1):
            if j != k:
                Augmented[k, j] = Augmented[k, j] / pivot_k
        
        # Actualizamos la columna pivote
        for i in range(n):
            if i != k:
                Augmented[i, k] = 0
        
        # El nuevo pivote es 1
        Augmented[k, k] = 1
        pivot = pivot_k

    # La solución está en la última columna
    x = Augmented[:, -1]
    return x

def Jacobi(A, b, tol=1e-10, max_iter=100):
    # Convertimos A y b a arrays de numpy
    A = np.array(A, dtype=float)
    b = np.array(b, dtype=float)
    
    # Dimensión del sistema
    n = len(b)
    
    # Vector inicializado en ceros
    x = np.zeros(n)
    
    # Iteración del método de Jacobi
    for k in range(max_iter):
        x_new = np.zeros_like(x)
        for i in range(n):
            suma = sum(A[i, j] * x[j] for j in range(n) if j != i)
            x_new[i] = (b[i] - suma) / A[i, i]
        
        # Verificamos la tolerancia
        if np.linalg.norm(x_new - x, ord=np.inf) < tol:
            return x_new
        
        x = x_new  # Actualizamos el vector para la siguiente iteración
    
    raise ValueError("El método de Jacobi no convergió en el número máximo de iteraciones.")

def GaussSeidel(A, b, tol=1e-10, max_iter=100):
    # Convertimos A y b a arrays de numpy
    A = np.array(A, dtype=float)
    b = np.array(b, dtype=float)
    
    # Dimensión del sistema
    n = len(b)
    
    # Vector inicializado en ceros
    x = np.zeros(n)
    
    # Iteración del método de Gauss-Seidel
    for k in range(max_iter):
        x_new = x.copy()
        for i in range(n):
            suma = sum(A[i, j] * x_new[j] for j in range(i)) + sum(A[i, j] * x[j] for j in range(i + 1, n))
            x_new[i] = (b[i] - suma) / A[i, i]
        
        # Verificamos la tolerancia
        if np.linalg.norm(x_new - x, ord=np.inf) < tol:
            return x_new
        
        x = x_new  # Actualizamos el vector para la siguiente iteración
    
    raise ValueError("El método de Gauss-Seidel no convergió en el número máximo de iteraciones.")

# Métodos de ecuaciones no lineales
# Método de Bisección
def biseccion(f, a, b, tol, max_iter):
    x = sp.symbols('x')
    f = sp.sympify(f)
    iter_count = 0
    while iter_count < max_iter:
        c = (a + b) / 2
        fc = f.subs(x, c)
        print(f"Iteración {iter_count + 1}: c = {c}, f(c) = {fc}")
        if abs(fc) < tol or abs(b - a) / 2 < tol:
            return c
        if f.subs(x, a) * fc < 0:
            b = c
        else:
            a = c
        iter_count += 1
    return c

# Método de Newton-Raphson
def newton_raphson(f, x0, tol, max_iter):
    x = sp.symbols('x')
    f = sp.sympify(f)
    df = sp.diff(f, x)  # Derivada de f
    iter_count = 0
    while iter_count < max_iter:
        fx = f.subs(x, x0)
        dfx = df.subs(x, x0)
        if dfx == 0:
            print("Derivada cero, no se puede continuar.")
            return None
        x1 = x0 - fx / dfx
        print(f"Iteración {iter_count + 1}: x = {x1}")
        if abs(x1 - x0) < tol:
            return x1
        x0 = x1
        iter_count += 1
    return x0

# Método de Falsa Posición
def falsa_posicion(f, a, b, tol, max_iter):
    x = sp.symbols('x')
    f = sp.sympify(f)
    iter_count = 0
    while iter_count < max_iter:
        fa = f.subs(x, a)
        fb = f.subs(x, b)
        c = b - (fb * (b - a)) / (fb - fa)
        fc = f.subs(x, c)
        print(f"Iteración {iter_count + 1}: c = {c}, f(c) = {fc}")
        if abs(fc) < tol or abs(b - a) < tol:
            return c
        if fa * fc < 0:
            b = c
        else:
            a = c
        iter_count += 1
    return c

# Método de Punto Fijo
def punto_fijo(g, x0, tol, max_iter):
    x = sp.symbols('x')
    g = sp.sympify(g)
    iter_count = 0
    while iter_count < max_iter:
        x1 = g.subs(x, x0)
        print(f"Iteración {iter_count + 1}: x = {x1}")
        if abs(x1 - x0) < tol:
            return x1
        x0 = x1
        iter_count += 1
    return x1

# Método de Secante
def secante(f, x0, x1, tol, max_iter):
    x = sp.symbols('x')
    f = sp.sympify(f)
    iter_count = 0
    while iter_count < max_iter:
        fx0 = f.subs(x, x0)
        fx1 = f.subs(x, x1)
        x2 = x1 - fx1 * (x1 - x0) / (fx1 - fx0)
        print(f"Iteración {iter_count + 1}: x = {x2}")
        if abs(x2 - x1) < tol:
            return x2
        x0, x1 = x1, x2
        iter_count += 1
    return x2

#Metodos de EDO
def euler(func, x0, y0, x_final, h):
    x, y = x0, y0
    resultados = []
    
    while x <= x_final:
        resultados.append(f"x = {x:.4f}, y = {y:.4f}")
        y += h * func(x, y)
        x += h

    return resultados  # Retorna los resultados formateados como cadenas


def runge_kutta_2(func, x0, y0, x_final, h):
    x, y = x0, y0
    resultados = []
    
    while x <= x_final:
        resultados.append(f"x = {x:.4f}, y = {y:.4f}")
        k1 = h * func(x, y)
        k2 = h * func(x + h, y + k1)
        y += (k1 + k2) / 2
        x += h

    return resultados  # Retorna los resultados formateados como cadenas


def runge_kutta_3(func, x0, y0, x_final, h):
    x, y = x0, y0
    resultados = []
    
    while x <= x_final:
        resultados.append(f"x = {x:.4f}, y = {y:.4f}")
        k1 = h * func(x, y)
        k2 = h * func(x + h / 2, y + k1 / 2)
        k3 = h * func(x + h, y - k1 + 2 * k2)
        y += (k1 + 4 * k2 + k3) / 6
        x += h

    return resultados  # Retorna los resultados formateados como cadenas

def euler_modificado(func, x0, y0, x_final, h):
    x, y = x0, y0
    resultados = []
    while x <= x_final:
        resultados.append(f"x = {x:.4f}, y = {y:.4f}")
        y_predict = y + h * func(x, y)
        y = y + h / 2 * (func(x, y) + func(x + h, y_predict))
        x += h
    return resultados  # Retorna los resultados formateados como cadenas

#Metodos de integracion
def regla_trapezoidal(funcion, a, b, n):
    """
    Calcula la integral utilizando la regla trapezoidal.
    """
    f = lambda x: eval(funcion)
    n = int(n)
    h = (b - a) / n
    suma = (f(a) + f(b)) / 2.0

    for i in range(1, n):
        x = a + i * h
        suma += f(x)

    integral = h * suma
    resultado = f"Resultado de la integral por la regla trapezoidal: {integral:.6f}"
    return resultado


def simpson_tercio(funcion, a, b, n):
    """
    Calcula la integral utilizando la regla de Simpson 1/3.
    """
    f = lambda x: eval(funcion)
    n = int(n)
    if n % 2 == 1:
        n += 1
    h = (b - a) / n
    suma = f(a) + f(b)

    for i in range(1, n, 2):
        suma += 4 * f(a + i * h)
    for i in range(2, n, 2):
        suma += 2 * f(a + i * h)

    integral = h / 3 * suma
    resultado = f"Resultado de la integral por la regla de Simpson 1/3: {integral:.6f}"
    return resultado


def simpson_octavos(funcion, a, b, n):
    """
    Calcula la integral utilizando la regla de Simpson 3/8.
    """
    f = lambda x: eval(funcion)
    n = int(n)
    if n % 3 != 0:
        n += 3 - (n % 3)
    h = (b - a) / n
    suma = f(a) + f(b)

    for i in range(1, n, 3):
        suma += 3 * f(a + i * h)
    for i in range(2, n, 3):
        suma += 3 * f(a + i * h)
    for i in range(3, n, 3):
        suma += 2 * f(a + i * h)

    integral = 3 * h / 8 * suma
    resultado = f"Resultado de la integral por la regla de Simpson 3/8: {integral:.6f}"
    return resultado


def newton_cotes_cerradas(funcion, a, b, n):
    """
    Calcula la integral utilizando Newton-Cotes cerradas.
    """
    f = lambda x: eval(funcion)
    n = int(n)
    h = (b - a) / n
    suma = 0

    for i in range(n + 1):
        coef = 1  # Placeholder for actual coefficients
        x = a + i * h
        suma += coef * f(x)

    integral = h * suma
    resultado = f"Resultado de la integral por Newton-Cotes cerradas: {integral:.6f}"
    return resultado


def newton_cotes_abiertas(funcion, a, b, n):
    """
    Calcula la integral utilizando Newton-Cotes abiertas.
    """
    f = lambda x: eval(funcion)
    n = int(n)
    h = (b - a) / (n + 1)
    suma = 0

    for i in range(1, n + 1):
        coef = 1  # Placeholder for actual coefficients
        x = a + i * h
        suma += coef * f(x)

    integral = h * suma
    resultado = f"Resultado de la integral por Newton-Cotes abiertas: {integral:.6f}"
    return resultado


#Minimos Cuadrados
# Línea Recta
def ajuste_linea_recta(x, y):
    n = len(x)
    x = np.array(x)
    y = np.array(y)
    m = (n * np.sum(x * y) - np.sum(x) * np.sum(y)) / (n * np.sum(x**2) - np.sum(x)**2)
    c = (np.sum(y) - m * np.sum(x)) / n
    return m, c

# Cuadrática
def ajuste_cuadratico(x, y):
    coeficientes = np.polyfit(x, y, 2)  # Ajuste polinómico de grado 2
    return coeficientes  # Devuelve [a, b, c] de ax^2 + bx + c

# Cúbica
def ajuste_cubico(x, y):
    coeficientes = np.polyfit(x, y, 3)  # Ajuste polinómico de grado 3
    return coeficientes  # Devuelve [a, b, c, d] de ax^3 + bx^2 + cx + d

# Lineal con Función (ejemplo con logaritmos)
def ajuste_lineal_funcion(x, y):
    x = np.array(x)
    y = np.log(np.array(y))  # Transformamos y con logaritmos
    m, c = ajuste_linea_recta(x, y)
    return m, np.exp(c)  # Convertimos c de nuevo con exponencial

# Cuadrática con Función (ejemplo con logaritmos)
def ajuste_cuadratico_funcion(x, y):
    x = np.array(x)
    y = np.log(np.array(y))  # Transformamos y con logaritmos
    coeficientes = np.polyfit(x, y, 2)  # Ajuste cuadrático en espacio transformado
    return coeficientes

# Función de cálculo de interpolación o ecuación
def calcular():
    try:
        if metodo_var.get() == "interpolacion":
            x = list(map(float, x_entry.get().split()))
            y = list(map(float, y_entry.get().split()))
            x_nuevo = float(x_nuevo_entry.get())
            
            if metodo_interpolacion_var.get() == "Lineal":
                resultado = interpolacion_lineal(x, y, x_nuevo)
            elif metodo_interpolacion_var.get() == "Lagrange":
                resultado = interpolacion_lagrange(x, y, x_nuevo)
            elif metodo_interpolacion_var.get() == "Newton Adelante":
                resultado = newton_adelante(x, y, x_nuevo)
            elif metodo_interpolacion_var.get() == "Newton Atras":
                resultado = newton_atras(x, y, x_nuevo)
            elif metodo_interpolacion_var.get() == "Newton con diferencias divididas":
                resultado = newton_diferencias_divididas(x, y, x_nuevo)                
            else:
                messagebox.showerror("Error", "Selecciona un método de interpolación válido.")
                
            messagebox.showinfo("Resultado", f"Resultado de la interpolación: {resultado}")
        
        elif metodo_var.get() == "lineal":
            A = np.array([list(map(float, row.split())) for row in A_entry.get().splitlines()])
            b = np.array(list(map(float, b_entry.get().split())))
            if metodo_lineal_var.get() == "Gauss-Jordan":
                resultado = gauss_jordan(A, b)
            elif metodo_lineal_var.get() == "Eliminación Gauss":
                resultado = gauss_eliminacion(A, b)
            elif metodo_lineal_var.get() == "Montante":
                resultado = Montante(A,b)
            elif metodo_lineal_var.get() == "Gauss-seidei":
                resultado = GaussSeidel(A,b)
            elif metodo_lineal_var.get() == "Jacobi":
                resultado = Jacobi(A,b)
            else:
                messagebox.showerror("Error", "Selecciona un método lineal válido.")
            
            messagebox.showinfo("Resultado", f"Solución del sistema: {resultado}")
        
        elif metodo_var.get() == "EDO":
            func_str = f_entry.get()
            func = lambda x, y: eval(func_str)
            x0 = float(x0_entry.get())
            y0 = float(y0_entry.get())
            x_final = float(x_final_entry.get())
            h = float(h_entry.get())
            if metodo_EDO_var.get() == "Euler Modificado":
                resultado = euler_modificado(func, x0, y0, x_final, h)
                for i, resultado in enumerate(resultado):
                    messagebox.showinfo("Resultado", f"Resultados de x{i + 1}, y{i + 1}:\n{resultado}")
            elif metodo_EDO_var.get() == "Runge-Kutta 2do Orden":
                resultado = runge_kutta_2(func, x0, y0, x_final, h)
                for i, resultado in enumerate(resultado):
                    messagebox.showinfo("Resultado", f"Resultados de x{i + 1}, y{i + 1}:\n{resultado}")
            elif metodo_EDO_var.get() == "Runge-Kutta 3er Orden":
                resultado = runge_kutta_3(func, x0, y0, x_final, h)
                for i, resultado in enumerate(resultado):
                    messagebox.showinfo("Resultado", f"Resultados de x{i + 1}, y{i + 1}:\n{resultado}")
            elif metodo_EDO_var.get() == "Euler":
                resultado =euler(func, x0, y0, x_final, h)
                for i, resultado in enumerate(resultado):
                    messagebox.showinfo("Resultado", f"Resultados de x{i + 1}, y{i + 1}:\n{resultado}")
            else:
              messagebox.showerror("Error", "Selecciona un método de EDO válido.")
              
            
        elif metodo_var.get() == "no_lineal":
            func_str = f_entry.get()  # La función se toma como input
            tol = float(tol_entry.get())  # Tolerancia
            max_iter = int(max_iter_entry.get())  # Número máximo de iteraciones
            f = lambda x: eval(func_str)  # Convertir la expresión en función

            if metodo_no_lineal_var.get() == "Bisección":
               a = float(a_entry.get())  # Límite inferior
               b = float(b_entry.get())  # Límite superior
               resultado = biseccion(f, a, b, tol, max_iter)
               messagebox.showinfo("Resultado", f"Raíz encontrada por Bisección: {resultado}")

            elif metodo_no_lineal_var.get() == "Newton-Raphson":
                x0 = float(x0_entry.get())  # Valor inicial
                resultado = newton_raphson(f, x0, tol, max_iter)
                messagebox.showinfo("Resultado", f"Raíz encontrada por Newton-Raphson: {resultado}")

            elif metodo_no_lineal_var.get() == "Falsa Posición":
                a = float(a_entry.get())  # Límite inferior
                b = float(b_entry.get())  # Límite superior
                resultado = falsa_posicion(f, a, b, tol, max_iter)
                messagebox.showinfo("Resultado", f"Raíz encontrada por Falsa Posición: {resultado}")

            elif metodo_no_lineal_var.get() == "Punto Fijo":
                g_str = g_entry.get()  # La función g(x) para el método de punto fijo
                g = lambda x: eval(g_str)
                x0 = float(x0_entry.get())  # Valor inicial
                resultado = punto_fijo(g, x0, tol, max_iter)
                messagebox.showinfo("Resultado", f"Raíz encontrada por Punto Fijo: {resultado}")

            elif metodo_no_lineal_var.get() == "Secante":
                x0 = float(x0_entry.get())  # Primer valor inicial
                x1 = float(x1_entry.get())  # Segundo valor inicial
                resultado = secante(f, x0, x1, tol, max_iter)
                messagebox.showinfo("Resultado", f"Raíz encontrada por Secante: {resultado}")

            else:
                messagebox.showerror("Error", "Selecciona un método no lineal válido.")
        
        elif metodo_var.get() == "Integracion":
            funcion = funcion_entry.get()  # La función se toma como input
            a = float(limite_inferior_entry.get())  # Límite inferior
            b = float(limite_superior_entry.get())  # Límite superior
            n = float(subintervalos_entry.get())  # Número de subintervalos

            if metodo_integracion_var.get() == "Regla Trapezoidal":
               resultado = regla_trapezoidal(funcion, a, b, n)
               messagebox.showinfo("Resultado", f"{resultado}")

            elif metodo_integracion_var.get() == "Regla de 1/3 Simpson":
                resultado = simpson_tercio(funcion, a, b, n)
                messagebox.showinfo("Resultado", f"{resultado}")

            elif metodo_integracion_var.get() == "Regla de 3/8 Simpson":
                resultado = simpson_octavos(funcion, a, b, n)
                messagebox.showinfo("Resultado", f"{resultado}")

            elif metodo_integracion_var.get() == "Newton-Cotes Cerradas":
                resultado = newton_cotes_cerradas(funcion, a, b, n)
                messagebox.showinfo("Resultado", f"{resultado}")

            elif metodo_integracion_var.get() == "Newton-Cotes Abiertas":
                resultado = newton_cotes_abiertas(funcion, a, b, n)
                messagebox.showinfo("Resultado", f"{resultado}")

            else:
                 messagebox.showerror("Error", "Selecciona un método válido de integración numérica.")
            return
        
        


        elif metodo_var.get() == "MinimosCuadrados":
            x = list(map(float, x_entry.get().split()))
            y = list(map(float, y_entry.get().split()))

            if ajuste_linea_recta == "Línea Recta":
                m, c = ajuste_linea_recta(x, y)
                resultado = f"Pendiente: {m}, Intersección: {c}"

            elif ajuste_cuadratico == "Cuadrática":
                a, b, c = ajuste_cuadratico(x, y)
                resultado = f"Coeficientes: a={a}, b={b}, c={c}"

            elif ajuste_cubico == "Cúbica":
                a, b, c, d = ajuste_cubico(x, y)
                resultado = f"Coeficientes: a={a}, b={b}, c={c}, d={d}"

            elif ajuste_lineal_funcion == "Lineal con Función":
                m, c = ajuste_lineal_funcion(x, y)
                resultado = f"Pendiente: {m}, Intersección (transformada): {c}"

            elif ajuste_cuadratico_funcion == "Cuadrática con Función":
                coeficientes = ajuste_cuadratico_funcion(x, y)
                resultado = f"Coeficientes (transformados): {coeficientes}"

            else:
                messagebox.showerror("Error", "Selecciona un método válido de mínimos cuadrados.")
                return
            

            messagebox.showinfo("Resultado", f"Resultado: {resultado}")
            

                    
    except Exception as e:
        messagebox.showerror("Error", f"Ocurrió un error: {str(e)}")


# Crear la interfaz de usuario (GUI)
root = tk.Tk()
root.title("Métodos Numéricos")

# Selección de tipo de problema
metodo_var = tk.StringVar(value="interpolacion")
tk.Label(root, text="Selecciona el tipo de problema:").pack()
metodo_menu = tk.OptionMenu(root, metodo_var, "interpolacion", "lineal", "no_lineal","EDO","Integracion","MinimosCuadrados", command=lambda x: cambiar_menu(x))
metodo_menu.pack()

# Contenido para interpolación
def cambiar_menu(tipo):
    # Ocultar todos los campos
    for widget in root.winfo_children():
        widget.pack_forget()

    # Menú principal
    tk.Label(root, text="Selecciona el tipo de problema:").pack()
    metodo_menu.pack()

    if tipo == "interpolacion":
        # Interpolación
        tk.Label(root, text="Selecciona el método de interpolación:").pack()
        interpolacion_menu.pack()
        tk.Label(root, text="Puntos x (separados por espacios):").pack()
        x_entry.pack()
        tk.Label(root, text="Puntos y (separados por espacios):").pack()
        y_entry.pack()
        tk.Label(root, text="Nuevo valor de x:").pack()
        x_nuevo_entry.pack()

    elif tipo == "lineal":
        # Ecuaciones Lineales
        tk.Label(root, text="Selecciona el método de ecuaciones lineales:").pack()
        lineal_menu.pack()
        tk.Label(root, text="Matriz A (filas separadas por saltos de línea):").pack()
        A_entry.pack()
        tk.Label(root, text="Vector b (separado por espacios):").pack()
        b_entry.pack()

    elif tipo == "no_lineal":
        # Ecuaciones No Lineales
        tk.Label(root, text="Selecciona el método de ecuaciones no lineales:").pack()
        no_lineal_menu.pack()
        tk.Label(root, text="Función f(x) (por ejemplo, x**2 - 2):").pack()
        f_entry.pack()
        tk.Label(root, text="Derivada f'(x) (para Newton-Raphson):").pack()
        df_entry.pack()
        tk.Label(root, text="Valor inicial x0:").pack()
        x0_entry.pack()
        tk.Label(root, text="Límite inferior (a) para Bisección/Falsa Posición:").pack()
        a_entry.pack()
        tk.Label(root, text="Límite superior (b) para Bisección/Falsa Posición:").pack()
        b_entry.pack()
        tk.Label(root, text="Tolerancia:").pack()
        tol_entry.pack()
        tk.Label(root, text="Número máximo de iteraciones:").pack()
        max_iter_entry.pack()

    elif tipo == "EDO":
        # Ecuaciones Diferenciales Ordinarias (EDO)
        tk.Label(root, text="Selecciona el método de ecuaciones diferenciales:").pack()
        EDO_menu.pack()
        tk.Label(root, text="Función f(x, y):").pack()
        f_entry.pack()
        tk.Label(root, text="Valor inicial x0:").pack()
        x0_entry.pack()
        tk.Label(root, text="Valor inicial y0:").pack()
        y0_entry.pack()
        tk.Label(root, text="Valor final x:").pack()
        x_final_entry.pack()
        tk.Label(root, text="Paso h:").pack()
        h_entry.pack()

    elif tipo == "Integracion":
        # Integración Numérica
        tk.Label(root, text="Selecciona el método de integración:").pack()
        metodo_integracion_menu.pack()
        tk.Label(root, text="Función a integrar (en términos de x):").pack()
        funcion_entry.pack()
        tk.Label(root, text="Límite inferior (a):").pack()
        limite_inferior_entry.pack()
        tk.Label(root, text="Límite superior (b):").pack()
        limite_superior_entry.pack()
        tk.Label(root, text="Número de subintervalos (n):").pack()
        subintervalos_entry.pack()

    elif tipo == "MinimosCuadrados":
        # Selección del tipo de ajuste
        tk.Label(root, text="Selecciona el tipo de ajuste:").pack()
        MinimosCuadrados_menu.pack()
        # Solicitar puntos x
        tk.Label(root, text="Puntos x (separados por espacios):").pack()
        x_minimos_entry = tk.Entry(root)
        x_minimos_entry.pack()
        # Solicitar puntos y
        tk.Label(root, text="Puntos y (separados por espacios):").pack()
        y_minimos_entry = tk.Entry(root)
        y_minimos_entry.pack()
    # Botón para calcular
    calc_button.pack()

# Opciones de interpolación
metodo_interpolacion_var = tk.StringVar(value="Lineal")
interpolacion_menu = tk.OptionMenu(root, metodo_interpolacion_var, "Lineal", "Lagrange", "Newton Atras", "Newton Adelante", "Newton con diferencias divididas")

# Opciones de ecuaciones lineales
metodo_lineal_var = tk.StringVar(value="Gauss-Jordan")
lineal_menu = tk.OptionMenu(root, metodo_lineal_var, "Gauss-Jordan", "Eliminación Gauss","Montante","Jacobi","Gauss-seidei")

# Opciones de ecuaciones no lineales
metodo_no_lineal_var = tk.StringVar(value="Bisección")
no_lineal_menu = tk.OptionMenu(root, metodo_no_lineal_var, "Bisectriz", "Newton-Raphson","Falsa-Posicion","Punto-Fijo","Secante")

# Opciones de EDO
metodo_EDO_var = tk.StringVar(value="Euler Modificado")
EDO_menu = tk.OptionMenu(root, metodo_EDO_var,
                         "Euler Modificado",
                         "Runge-Kutta 2do Orden",
                         "Runge-Kutta 3er Orden",
                         "Euler")

#Opciones Integracion
metodo_integracion_var = tk.StringVar(value="Regla Trapezoidal")  # Valor por defecto
metodo_integracion_menu = tk.OptionMenu(root,metodo_integracion_var,
    "Regla Trapezoidal",
    "Regla de 1/3 Simpson",
    "Regla de 3/8 Simpson",
    "Newton-Cotes Cerradas",
    "Newton-Cotes Abiertas"
    )

# Opciones MinimosCuadrados
metodo_MinimosCuadrados_var = tk.StringVar(value="Línea Recta")  # Valor por defecto
MinimosCuadrados_menu = tk.OptionMenu(root,metodo_MinimosCuadrados_var,
    "Línea Recta",             # Ajuste lineal
    "Cuadrática",              # Ajuste polinómico de grado 2
    "Cúbica",                  # Ajuste polinómico de grado 3
    "Lineal con Función",      # Ajuste lineal transformado (ejemplo: logaritmo)
    "Cuadrática con Función"   # Ajuste cuadrático transformado (ejemplo: logaritmo)
)

# Entradas y botón de cálculo
x_entry = tk.Entry(root)
y_entry = tk.Entry(root)
x_nuevo_entry = tk.Entry(root)
A_entry = tk.Entry(root)
b_entry = tk.Entry(root)
f_entry = tk.Entry(root)
df_entry = tk.Entry(root)
x0_entry = tk.Entry(root)
x0_entry = tk.Entry(root)
y0_entry = tk.Entry(root)
x_final_entry = tk.Entry(root)
h_entry = tk.Entry(root)
funcion_entry = tk.Entry(root)
limite_inferior_entry = tk.Entry(root)
limite_superior_entry = tk.Entry(root)
subintervalos_entry = tk.Entry(root)
a_entry = tk.Entry(root)
tol_entry = tk.Entry(root)
max_iter_entry = tk.Entry(root)
x1_entry = 0
g_entry = 0
calc_button = tk.Button(root, text="Calcular", command=calcular)

# Iniciar la interfaz gráfica
root.mainloop()
