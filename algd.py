a = "ATTCGTCCAGTCCTGA"
b = "TTACGGCGTAACCTGGTA"
NB = []
PB = []

def arreglosimbolos(a,b):
    simbolos = []
    lista = list(a)
    listb = list(b)
    for i in lista:
        valor = 0
        for j in simbolos:
            if (i == j):
                valor = 1
        if (valor == 0):
            simbolos.append(i)
      
    
    for i in listb:
        valor = 0
        for j in simbolos:
            if (i == j):
                valor = 1
        if (valor == 0):
            simbolos.append(i)
                
    print simbolos
    return simbolos
                
                
def cuantasimbolopalabra(palabra,simbolo):
    lista = list(palabra)
    contador = 0
    for i in lista:
        if (simbolo == i):
            contador = contador + 1
            
            
    print contador
    return contador



def posicionessimbolo(palabra,simbolo,contador):
    lista = list(palabra)
    pb = []
    for i in range(0,len(lista)):
        
        if (lista[i] == simbolo):
            
            pb.append(i)

    print pb
    return pb

simbolos = arreglosimbolos(a,b)
contador = cuantasimbolopalabra(a,"G")
posiciones = posicionessimbolo(a,"T",contador)