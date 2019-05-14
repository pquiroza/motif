import matplotlib.pyplot as plt

def leearchivo(nombre):
    with open(nombre,"r") as ins:
        array = []
        for line in ins:
            array.append(float(line))
    
    import matplotlib.pyplot as plt
    plt.bar(range(len(array)),array)
    plt.show()
    
leearchivo("datos.fts")