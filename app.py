from flask import Flask, render_template, request, jsonify
import numpy as np
from sympy import symbols, sympify, lambdify, diff
import re

app = Flask(__name__)

def convert_power_operator(function):
    function = function.replace('^', '**')
    function = function.replace(')(', ')*(')
    function = re.sub(r'(\d)([a-zA-Z])', r'\1*\2', function)
    function = re.sub(r'([a-zA-Z])(\d)', r'\1*\2', function)
    function = re.sub(r'([a-zA-Z0-9])\(', r'\1*(', function)
    function = re.sub(r'\)([a-zA-Z0-9])', r')*\1', function)
    return function

def create_function(f_str):
    x = symbols('x')
    f_str = convert_power_operator(f_str)
    try:
        f_expr = sympify(f_str)
    except Exception as e:
        raise ValueError(f"Format fungsi tidak valid: {str(e)}")

    f_lambda = lambdify(x, f_expr, modules=['numpy'])
    
    def f_numeric(x_val):
        try:
            result = float(f_lambda(x_val))
            if not np.isfinite(result):
                raise ValueError("Hasil perhitungan tidak valid")
            return result
        except Exception as e:
            raise ValueError(f"Error dalam perhitungan f({x_val}): {str(e)}")
    
    return f_expr, f_numeric

def interval_halving(f, x0, x1, epsilon):
    iterations = []
    try:
        f_expr, f_numeric = create_function(f)
        
        a = float(x0)  # Batas bawah
        b = float(x1)  # Batas atas
        
        if a >= b:
            raise ValueError("Interval awal (a) harus lebih kecil dari interval akhir (b)")

        fa = f_numeric(a)   # f(a)
        fb = f_numeric(b)   # f(b)
        
        if fa * fb >= 0:
            raise ValueError("Interval awal tidak memenuhi syarat teorema (f(a)*f(b) harus < 0)")
        
        iteration = 0
        while True:
            c = (a + b) / 2  # Titik tengah
            fa = f_numeric(a)   # f(a)
            fb = f_numeric(b)   # f(b)
            fc = f_numeric(c)   # f(c)
            
            iteration_data = {
                'iteration': iteration + 1,
                'a': a,
                'b': b,
                'c': c,
                'f(a)': fa,
                'f(b)': fb,
                'f(c)': fc,
                'epsilon': epsilon
            }
            iterations.append(iteration_data)
            
            if abs(fc) < epsilon:
                break
                
            if fa * fc < 0:
                b = c
            else:
                a = c
                
            iteration += 1
            if iteration > 100:
                raise ValueError("Iterasi maksimum tercapai (100 iterasi)")
        
        return iterations, c
        
    except Exception as e:
        raise ValueError(f"Error dalam perhitungan: {str(e)}")

def regula_falsi(f, x0, x1, epsilon):
    iterations = []
    try:
        f_expr, f_numeric = create_function(f)
        
        xn = float(x0)
        xn1 = float(x1)
        
        if xn >= xn1:
            raise ValueError("Interval awal (xₙ) harus lebih kecil dari interval akhir (xₙ₊₁)")

        fxn = f_numeric(xn)
        fxn1 = f_numeric(xn1)
        
        if fxn * fxn1 >= 0:
            raise ValueError("Interval awal tidak memenuhi syarat teorema (f(xₙ)*f(xₙ₊₁) harus < 0)")
        
        iteration = 0
        while True:
            # Rumus Regula Falsi sesuai flowchart
            xt = xn - fxn * (xn1 - xn) / (fxn1 - fxn)
            fxt = f_numeric(xt)
            
            iteration_data = {
                'iteration': iteration + 1,
                'xn': xn,
                'xn1': xn1,
                'xt': xt,
                'f(xn)': fxn,
                'f(xn1)': fxn1,
                'f(xt)': fxt,
                'epsilon': epsilon
            }
            iterations.append(iteration_data)
            
            if abs(fxt) < epsilon:
                break
                
            if fxn * fxt < 0:
                xn1 = xt
                fxn1 = fxt
            else:
                xn = xt
                fxn = fxt
                
            iteration += 1
            if iteration > 100:
                raise ValueError("Iterasi maksimum tercapai (100 iterasi)")
        
        return iterations, xt
        
    except Exception as e:
        raise ValueError(f"Error dalam perhitungan: {str(e)}")

def newton_raphson(f, x0, epsilon):
    iterations = []
    try:
        x = symbols('x')
        f_expr, f_numeric = create_function(f)
        
        # Langkah 1: Inisialisasi
        x1 = float(x0)
        iteration = 0
        
        # Langkah 2: Hitung turunan
        df_expr = diff(f_expr, x)  # Ekspresi turunan
        df_lambda = lambdify(x, df_expr)
        
        # Langkah 3: Iterasi
        while True:
            fx1 = f_numeric(x1)
            dfx = df_lambda(x1)
            
            if abs(dfx) < 1e-10:
                raise ValueError("f'(x) terlalu kecil")
            
            x2 = x1 - fx1/dfx  # Rumus Newton
            fx2 = f_numeric(x2)
            
            iteration_data = {
                'iteration': iteration + 1,
                'x1': x1,
                'f(x1)': fx1,
                'f\'(x)': dfx,
                'df_expr': str(df_expr),  # Tambahkan ekspresi turunan
                'x2': x2,
                'f(x2)': fx2,
                'epsilon': epsilon
            }
            iterations.append(iteration_data)
            
            if abs(fx2) < epsilon:  # Kriteria berhenti
                break
                
            x1 = x2
            iteration += 1
            if iteration > 100:
                raise ValueError("Iterasi maksimum tercapai (100 iterasi)")
        
        return iterations, x2
        
    except Exception as e:
        raise ValueError(f"Error dalam perhitungan: {str(e)}")

def secant(f, x0, x1, epsilon):
    iterations = []
    try:
        f_expr, f_numeric = create_function(f)
        
        x_prev = float(x0)
        x_curr = float(x1)
        
        # Evaluasi fungsi di titik awal
        f_prev = f_numeric(x_prev)
        f_curr = f_numeric(x_curr)
        
        iteration = 0
        while True:
            # Cek pembagian dengan nol
            if abs(f_curr - f_prev) < 1e-10:
                raise ValueError("Pembagian dengan nilai terlalu kecil")
            
            # Hitung x berikutnya menggunakan rumus secant
            x_next = x_curr - f_curr * (x_curr - x_prev) / (f_curr - f_prev)
            f_next = f_numeric(x_next)
            
            iteration_data = {
                'iteration': iteration + 1,
                'x0': x_prev,
                'x1': x_curr,
                'f(x0)': f_prev,
                'f(x1)': f_curr,
                'x2': x_next,
                'f(x2)': f_next,
                'epsilon': epsilon
            }
            iterations.append(iteration_data)
            
            # Kriteria konvergensi: selisih antara dua titik berurutan
            if abs(x_next - x_curr) < epsilon:
                break
            
            # Update nilai untuk iterasi berikutnya
            x_prev = x_curr
            x_curr = x_next
            f_prev = f_curr
            f_curr = f_next
            
            iteration += 1
            if iteration > 100:
                raise ValueError("Iterasi maksimum tercapai (100 iterasi)")
        
        return iterations, x_next
        
    except Exception as e:
        raise ValueError(f"Error dalam perhitungan: {str(e)}")

@app.route('/')
def home():
    return render_template('index.html')

@app.route('/calculate', methods=['POST'])
def calculate():
    try:
        data = request.get_json()
        function = data['function']
        method = data['method']
        epsilon = float(data['epsilon'])
        
        if method == 'bisection':
            x0 = float(data['x0'])
            x1 = float(data['x1'])
            iterations, root = interval_halving(function, x0, x1, epsilon)
        elif method == 'regulafalsi':
            x0 = float(data['x0'])
            x1 = float(data['x1'])
            iterations, root = regula_falsi(function, x0, x1, epsilon)
        elif method == 'newton':
            x0 = float(data['x0'])
            iterations, root = newton_raphson(function, x0, epsilon)
        elif method == 'secant':
            x0 = float(data['x0'])
            x1 = float(data['x1'])
            iterations, root = secant(function, x0, x1, epsilon)
        else:
            raise ValueError("Metode tidak valid")
        
        return jsonify({
            'success': True,
            'iterations': iterations,
            'root': root
        })
    except Exception as e:
        return jsonify({
            'success': False,
            'error': str(e)
        })

if __name__ == '__main__':
    app.run(host='0.0.0.0', port=5000, debug=True) 