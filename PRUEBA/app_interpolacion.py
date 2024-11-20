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
        else:
            messagebox.showerror("Error", "Selecciona un tipo de problema válido.")
        
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
    
    elif tipo =="EDO":
        tk.Label(root, text="Selecciona el método de ecuaciones diferenciales:").pack()
        EDO_menu.pack()

    elif tipo =="Integracion":
        tk.Label(root, text="Selecciona el método de integracion:").pack()
        Integracion_menu.pack()

    elif tipo =="Minimos cuadrados":
        tk.Label(root, text="Selecciona el método de minimos cuadrados:").pack()
        MinimosCuadrados_menu.pack()

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
metodo_EDO_var = tk.StringVar(value="Euler_Modificado")
EDO_menu = tk.OptionMenu(root, metodo_EDO_var,"EulerModificado")

#Opciones Integracion
metodo_Integracion_var = tk.StringVar(value="N-Cotes")
Integracion_menu = tk.OptionMenu(root, metodo_Integracion_var,"N-Cotes")

#Opciones MinimosCuadrados
metodo_MinimosCuadrados_var = tk.StringVar(value="Linea-Recta")
MinimosCuadrados_menu = tk.OptionMenu(root, metodo_MinimosCuadrados_var,"Linea-recta")

# Entradas y botón de cálculo
x_entry = tk.Entry(root)
y_entry = tk.Entry(root)
x_nuevo_entry = tk.Entry(root)
A_entry = tk.Entry(root)
b_entry = tk.Entry(root)
f_entry = tk.Entry(root)
df_entry = tk.Entry(root)
x0_entry = tk.Entry(root)
calc_button = tk.Button(root, text="Calcular", command=calcular)

# Iniciar la interfaz gráfica
root.mainloop()
