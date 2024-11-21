import tkinter as tk
from tkinter import messagebox
import numpy as np
from scipy.interpolate import lagrange

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

# Newton Adelante
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

# Métodos de ecuaciones no lineales
def biseccion(f, a, b, tol=1e-5, max_iter=100):
    if f(a) * f(b) > 0:
        raise ValueError("f(a) y f(b) deben tener signos opuestos")
    for _ in range(max_iter):
        c = (a + b) / 2
        if abs(f(c)) < tol:
            return c
        elif f(c) * f(a) < 0:
            b = c
        else:
            a = c
    return c

def newton_raphson(f, df, x0, tol=1e-5, max_iter=100):
    x = x0
    for _ in range(max_iter):
        fx = f(x)
        dfx = df(x)
        if abs(fx) < tol:
            return x
        x = x - fx / dfx
    return x

#Metodos de EDO
def runge_kutta_4(func, x0, y0, x_final, h):
    x, y = x0, y0
    resultados = []
    while x <= x_final:
        resultados.append((x, y))
        k1 = h * func(x, y)
        k2 = h * func(x + h / 2, y + k1 / 2)
        k3 = h * func(x + h / 2, y + k2 / 2)
        k4 = h * func(x + h, y + k3)
        y += (k1 + 2 * k2 + 2 * k3 + k4) / 6
        x += h
    return resultados

def runge_kutta_3_8(func, x0, y0, x_final, h):
    x, y = x0, y0
    resultados = []
    while x <= x_final:
        resultados.append((x, y))
        k1 = h * func(x, y)
        k2 = h * func(x + h / 3, y + k1 / 3)
        k3 = h * func(x + 2 * h / 3, y - k1 / 3 + k2)
        k4 = h * func(x + h, y + k1 - k2 + k3)
        y += (k1 + 3 * k2 + 3 * k3 + k4) / 8
        x += h
    return resultados

def runge_kutta_1_3(func, x0, y0, x_final, h):
    x, y = x0, y0
    resultados = []
    while x <= x_final:
        resultados.append((x, y))
        k1 = h * func(x, y)
        k2 = h * func(x + h / 2, y + k1 / 2)
        k3 = h * func(x + h, y - k1 + 2 * k2)
        y += (k1 + 4 * k2 + k3) / 6
        x += h
    return resultados

def runge_kutta_2(func, x0, y0, x_final, h):
    x, y = x0, y0
    resultados = []
    while x <= x_final:
        resultados.append((x, y))
        k1 = h * func(x, y)
        k2 = h * func(x + h, y + k1)
        y += (k1 + k2) / 2
        x += h
    return resultados

def runge_kutta_3(func, x0, y0, x_final, h):
    x, y = x0, y0
    resultados = []
    while x <= x_final:
        resultados.append((x, y))
        k1 = h * func(x, y)
        k2 = h * func(x + h / 2, y + k1 / 2)
        k3 = h * func(x + h, y - k1 + 2 * k2)
        y += (k1 + 4 * k2 + k3) / 6
        x += h
    return resultados

def euler_modificado(func, x0, y0, x_final, h):
    x, y = x0, y0
    resultados = []
    while x <= x_final:
        resultados.append((x, y))
        y_predict = y + h * func(x, y)
        y = y + h / 2 * (func(x, y) + func(x + h, y_predict))
        x += h
    return resultados

#Metodos de integracion
def regla_trapezoidal(funcion, a, b, n):
    # Convertir la función en un objeto de Python (usando eval, para simplificación)
    f = lambda x: eval(funcion)

    # Calcular el tamaño de cada subintervalo
    h = (b - a) / n

    # Aproximación inicial (suma de los extremos)
    suma = (f(a) + f(b)) / 2.0

    # Sumar las evaluaciones intermedias
    for i in range(1, n):
        x = a + i * h
        suma += f(x)

    # Multiplicar por el tamaño del intervalo
    integral = h * suma
    return integral

def simpson_tercio(funcion, a, b, n):
    # Convertir la función en un objeto de Python
    f = lambda x: eval(funcion)

    # Verificar que el número de subintervalos sea par
    if n % 2 == 1:
        n += 1  # Si es impar, hacer el número de subintervalos par

    # Calcular el tamaño de cada subintervalo
    h = (b - a) / n

    # Aproximación inicial
    suma = f(a) + f(b)

    # Sumar las evaluaciones de los puntos impares y pares
    for i in range(1, n, 2):
        suma += 4 * f(a + i * h)
    for i in range(2, n-1, 2):
        suma += 2 * f(a + i * h)

    # Multiplicar por h/3
    integral = h / 3 * suma
    return integral

def simpson_octavos(funcion, a, b, n):
    # Convertir la función en un objeto de Python
    f = lambda x: eval(funcion)

    # Verificar que el número de subintervalos sea múltiplo de 3
    if n % 3 != 0:
        n += 3 - (n % 3)  # Ajustar para que sea múltiplo de 3

    # Calcular el tamaño de cada subintervalo
    h = (b - a) / n

    # Aproximación inicial
    suma = f(a) + f(b)

    # Sumar las evaluaciones de los puntos
    for i in range(1, n, 3):
        suma += 3 * f(a + i * h)
    for i in range(2, n-1, 3):
        suma += 3 * f(a + i * h)
    for i in range(3, n-1, 3):
        suma += 2 * f(a + i * h)

    # Multiplicar por 3h/8
    integral = 3 * h / 8 * suma
    return integral

def newton_cotes_cerradas(funcion, a, b, n):
    # Convertir la función en un objeto de Python
    f = lambda x: eval(funcion)

    # Usar la fórmula de Newton-Cotes cerrada (generalizada)
    # Se necesita una fórmula específica para el valor de n (n puntos)

    h = (b - a) / n
    suma = 0

    for i in range(n + 1):
        coef = 1  # Esto cambiaría dependiendo de los coeficientes para cada n (simplificado aquí)
        x = a + i * h
        suma += coef * f(x)

    integral = h * suma  # Ajustar según los coeficientes de la fórmula
    return integral

def newton_cotes_abiertas(funcion, a, b, n):
    # Convertir la función en un objeto de Python
    f = lambda x: eval(funcion)

    # Aquí necesitaríamos un procedimiento específico para los métodos de Newton-Cotes abiertos
    # Se utiliza una fórmula general de Newton-Cotes para puntos internos del intervalo.

    h = (b - a) / (n + 1)  # Ajustado para los puntos internos
    suma = 0

    for i in range(1, n + 1):
        coef = 1  # Esto también cambiaría dependiendo del número de puntos (simplificado aquí)
        x = a + i * h
        suma += coef * f(x)

    integral = h * suma  # Ajuste según los coeficientes
    return integral

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
                resultado = lagrange(x, y, x_nuevo)
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
            elif metodo_EDO_var.get() == "Runge-Kutta 2do Orden":
                resultado = runge_kutta_2(func, x0, y0, x_final, h)
            elif metodo_EDO_var.get() == "Runge-Kutta 3er Orden":
                resultado = runge_kutta_3(func, x0, y0, x_final, h)
            elif metodo_EDO_var.get() == "Runge-Kutta 4to Orden":
                resultado = runge_kutta_4(func, x0, y0, x_final, h)
            elif metodo_EDO_var.get() == "Runge-Kutta 3/8 Simpson":
                resultado = runge_kutta_3_8(func, x0, y0, x_final, h)
            elif metodo_EDO_var.get() == "Runge-Kutta 1/3 Simpson":
                resultado = runge_kutta_1_3(func, x0, y0, x_final, h)
            else:
              messagebox.showerror("Error", "Selecciona un método de EDO válido.")
              
              messagebox.showinfo("Resultado", f"Resultados: {resultado}")
            
        elif metodo_var.get() == "no_lineal":
            func_str = f_entry.get()
            f = lambda x: eval(func_str)  # Convertir la expresión en función
            df_str = df_entry.get()
            df = lambda x: eval(df_str)
            x0 = float(x0_entry.get())
            
            if metodo_no_lineal_var.get() == "Bisectriz":
                resultado = biseccion(f, -10, 10)  # Intervalo de ejemplo
            elif metodo_no_lineal_var.get() == "Newton-Raphson":
                resultado = newton_raphson(f, df, x0)
            else:
                messagebox.showerror("Error", "Selecciona un método no lineal válido.")
            
            messagebox.showinfo("Resultado", f"Raíz encontrada: {resultado}")

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
        
        elif n <= 0:
            funcion = funcion_entry.get()
            a = float(limite_inferior_entry.get())
            b = float(limite_superior_entry.get())
            n = int(subintervalos_entry.get())
            metodo = metodo_integracion_var.get()
            raise ValueError("El número de subintervalos debe ser mayor que cero.")

            # Realizar cálculo según el método seleccionado
            if metodo == "Regla Trapezoidal":
                resultado = regla_trapezoidal(funcion, a, b, n)
            elif metodo == "Regla de 1/3 Simpson":
                resultado = simpson_tercio(funcion, a, b, n)
            elif metodo == "Regla de 3/8 Simpson":
                resultado = simpson_octavos(funcion, a, b, n)
            elif metodo == "Newton – Cotes Cerradas":
                resultado = newton_cotes_cerradas(funcion, a, b, n)
            elif metodo == "Newton – Cotes Abiertas":
                resultado = newton_cotes_abiertas(funcion, a, b, n)
            else:
                resultado = "Método no reconocido."

            messagebox.showinfo("Resultado", f"Resultado de la integración: {resultado}")
                    
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

    elif tipo == "EDO":
        tk.Label(root, text="Selecciona el método de ecuaciones diferenciales:").pack()
        EDO_menu.pack()

        # Campos para ingresar la función, valores iniciales y parámetros
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


    elif tipo == "Integración":
        # Selección del método de integración
        tk.Label(root, text="Selecciona el método de integración:").pack()
        # Solicitar función a integrar
        tk.Label(root, text="Función a integrar (en términos de x):").pack()
        funcion_entry = tk.Entry(root)
        funcion_entry.pack()

        # Solicitar límites de integración
        tk.Label(root, text="Límite inferior (a):").pack()
        limite_inferior_entry = tk.Entry(root)
        limite_inferior_entry.pack()

        tk.Label(root, text="Límite superior (b):").pack()
        limite_superior_entry = tk.Entry(root)
        limite_superior_entry.pack()

        # Solicitar número de subintervalos
        tk.Label(root, text="Número de subintervalos (n):").pack()
        subintervalos_entry = tk.Entry(root)
        subintervalos_entry.pack()



    elif tipo == "Minimos cuadrados":
        # Selección del tipo de ajuste
        tk.Label(root, text="Selecciona el tipo de ajuste:").pack()
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
                         "Runge-Kutta 4to Orden",
                         "Runge-Kutta 3/8 Simpson",
                         "Runge-Kutta 1/3 Simpson")

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
                         "Runge-Kutta 4to Orden",
                         "Runge-Kutta 3/8 Simpson",
                         "Runge-Kutta 1/3 Simpson")

#Opciones Integracion
metodo_integracion_var = tk.StringVar(value="Regla Trapezoidal")  # Valor por defecto
metodo_integracion_menu = tk.OptionMenu(
    root,
    metodo_integracion_var,
    "Regla Trapezoidal",
    "Regla de 1/3 Simpson",
    "Regla de 3/8 Simpson",
    "Newton – Cotes Cerradas",
    "Newton – Cotes Abiertas"
    )

# Opciones MinimosCuadrados
metodo_MinimosCuadrados_var = tk.StringVar(value="Línea Recta")  # Valor por defecto

# Menú desplegable con todas las opciones
MinimosCuadrados_menu = tk.OptionMenu(
    root,
    metodo_MinimosCuadrados_var,
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
limite_inferior_entry = tk.Entry
limite_superior_entry = tk.Entry
subintervalos_entry = tk.Entry
calc_button = tk.Button(root, text="Calcular", command=calcular)

# Iniciar la interfaz gráfica
root.mainloop()
