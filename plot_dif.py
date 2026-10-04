import pandas as pd
import matplotlib.pyplot as plt
import os



files = [
    "func1_sin_x2_errors.csv",
    "func2_cos_sin_x_errors.csv",
    "func3_exp_sin_cos_x_errors.csv",
    "func4_ln_x_plus_3_errors.csv",
    "func5_sqrt_x_plus_3_errors.csv"
]

titles = [
    r"$f(x) = \sin(x^2)$",
    r"$f(x) = \cos(\sin(x))$",
    r"$f(x) = \exp(\sin(\cos(x)))$",
    r"$f(x) = \ln(x+3)$",
    r"$f(x) = (x+3)^{0.5}$"
]



method_labels = [
    "Правая разность",
    "Левая разность",
    "Центральная разность",
    "4й порядок",
    "6й порядок"
]

colors = ['blue', 'orange', 'green', 'red', 'purple']

for file, title in zip(files, titles):

    df = pd.read_csv(file)

    plt.figure(figsize=(10, 6))

    for i, col in enumerate(['method1', 'method2', 'method3', 'method4', 'method5']):
        plt.loglog(df['h'], df[col], 
                   marker="o", 
                   linestyle='-', 
                   color=colors[i], 
                   label=method_labels[i],
                   markersize=4)

    plt.xlabel('Шаг', fontsize=12)
    plt.ylabel('Абсолютная погрешность', fontsize=12)
    plt.title(f'{title}', fontsize=14)
    
    plt.grid(True, which="both", ls="--", alpha=0.6)
    plt.legend(fontsize=10)
    
    plt.tight_layout()
    plt.show()
