import random
import copy
import decimal
import sys
import uuid
import hashlib
import time
import math


archivo = sys.argv[3]
poblacion = 500000
largo = int(sys.argv[1])
datos = []
pop = []

class Kmer:
    def __init__(self,idgenome,adn  ,posicioninicial,largo):
        self.idgenome=idgenome
        self.adn=adn
        self.posicioninicial = posicioninicial
        self.largo = largo

class nmatrix:
    def __init__(self, indices,fitness,code):
        self.indices = indices
        self.fitness = fitness
        self.code = code

    def setindice(self,indices):
        self.indices = indices
    def setFitness(self,fitness):
        self.fitness = fitness
    def setCode(self,code):
        self.code = code

class matrix:
    def __init__(self, kmers,fitness,fitness2,motif,code):
        self.kmers = kmers
        self.fitness = fitness
        self.fitness2 = fitness2
        self.motif = motif
        self.code = code

    def setCode(self,code):
        self.code = code

    def setKmer(self,posicion,kmer):
        self.kmers[posicion] = kmer

    def setFitness(self,fitness):
        self.fitness = fitness

    def setFitness2(self,fitness2):
        self.fitness2 = fitness2

    def setMotif(self,motif):
        self.motif = motif


def showMatriz(matriz):
    for i in matriz.kmers:
        print (i.adn)
    print (matriz.fitness)

def escribeFitness(fitness,lkmer,archivo):
    f = open(archivo,'a')
    f.write(str(fitness))
    f.write('\n')
    f.close



def guardaGeneracion(archivo,pop):

    f = open(archivo,'w')
    for p in pop:
        pal = ""
        for i in p.indices:
            pal = pal+str(i)+","
        f.write(pal)
        f.write('\n')
    f.close

def escribeindividuo(individuo,generacion,archivo):
    f=open(archivo,'w')
    f.write("Motifs "+str(generacion))
    f.write('\n')
    ind = 0
    for i in individuo.indices:

        f.write(datos[ind][i:i+largo] +" " +str(i))
        f.write('\n')
        ind = ind + 1


    f.write(str(individuo.fitness))
    f.write('\n')
    f.close

def cargaRespaldo(archivo):
    print (archivo)
    poprespaldada = []
    f = open(archivo,'r')
    for line in f:
        d = line.split(",")
        results = d[:len(d)-1]
        results = [int(i) for i in results]
        m = nmatrix(results,0,"")
        poprespaldada.append(m)
    return poprespaldada

def escribeindividuo(individuo,generacion,archivo):
    f=open(archivo,'w')
    f.write("Motifs "+str(generacion))
    f.write('\n')
    ind = 0
    for i in individuo.indices:

        f.write(datos[ind][i:i+largo] +" " +str(i))
        f.write('\n')
        ind = ind + 1


    f.write(str(individuo.fitness))
    f.write('\n')
    f.close

def cargaIndividuo(archivo):
    datos = []
    f = open(archivo,'r')
    for line in f:
        valor = int(line)
        datos.append(valor)

    matriz = nmatrix(datos,0,"")
    return matriz

def cargaSecuencias(archivo):

    datos = []
    f = open(archivo,'r')
    secuencia = ""
    for line in f:
        p = line
        if p[0]==">":
                if (secuencia!=""):
                    datos.append(secuencia)
                    secuencia = ""
                else:
                    sp = p.split("|",2)
                    ns = sp[2].split(" ",2)
                    nombre = ns[0]

                    secuencia = ""

        else:
            l = line.split("\n",1)

            secuencia = secuencia + l[0]
    datos.append(secuencia)



    return datos





def generaMatrizInicialNuevo(largo):

    kmers = []
    matriz = nmatrix([],0,"")
    for i in datos:
        inicio = random.randint(0,len(i)-largo)
        kmers.append(inicio)
    matriz.setindice(kmers)
    matriz.setFitness(0)
    return matriz


def AnnealingFitness(matriz):
    largodatos = len(datos)
    alfabeto=["A","C","E","D","G","F","I","H","K","M","L","N","Q","P","S","R","T","W","V","Y","X"]
    ind = 0
    palabras = []
    fitness = 0
    for k in matriz.indices:
        palabras.append(datos[ind][k:k+largo])
        ind = ind + 1
    contador = 1
    for i in range(largo):
        counts = {"A":0,"C":0,"D":0,"E":0,"F":0,"G":0,"H":0,"I":0,"K":0,"L":0,"M":0,"N":0,"P":0,"Q":0,"R":0,"S":0,"T":0,"V":0,"W":0,"Y":0,"X":0}
        for l in palabras:
            counts[l[i]]+=1
        valor = counts.values()
        maximo = max(valor)
        if(maximo == largodatos):

            maximo = maximo + largo

        fitness = fitness + maximo
        ceros = list(counts.values()).count(0)
        fitness = fitness + ceros

    propor = decimal.Decimal(decimal.Decimal(fitness)/(decimal.Decimal(len(datos)+20+largo)*largo))

    matriz.setFitness(propor)



def newFitness(mmatriz):
    localr = []
    largodatos = len(datos)
    alfabeto=["A","C","E","D","G","F","I","H","K","M","L","N","Q","P","S","R","T","W","V","Y","X"]
    for matriz in mmatriz:
        ind = 0
        palabras = []
        fitness=0
        for k in matriz.indices:

            palabras.append(datos[ind][k:k+largo])
            ind = ind + 1
        contador = 1
        for i in range(largo):

                counts = {"A":0,"C":0,"D":0,"E":0,"F":0,"G":0,"H":0,"I":0,"K":0,"L":0,"M":0,"N":0,"P":0,"Q":0,"R":0,"S":0,"T":0,"V":0,"W":0,"Y":0,"X":0}

                for l in palabras:
                    counts[l[i]]+=1
                #print(counts)

                valor = counts.values()
                maximo = max(valor)

                if(maximo == largodatos):

                    maximo = maximo + largo

                fitness = fitness + maximo
                ceros = list(counts.values()).count(0)
                fitness = fitness + ceros



        propor = decimal.Decimal(decimal.Decimal(fitness)/(decimal.Decimal(len(datos)+20+largo)*largo))

        matriz.setFitness(propor)



def mutacion2(m):
    vecino = copy.deepcopy(m)
    for i in range(len(datos)):
        gen = random.randint(0,len(vecino.indices)-1)
        nuevogen = random.randint(0,len(datos[gen])-largo)
        vecino.indices[gen] = nuevogen
    return vecino

def mutacion(m):
    vecino = copy.deepcopy(m)
    gen = random.randint(0,len(vecino.indices)-1)
    nuevogen = random.randint(0,len(datos[gen])-largo)
    vecino.indices[gen] = nuevogen

    return vecino

def cruzamiento(m1,m2):

    hijo1 = copy.deepcopy(m1)
    hijo2 = copy.deepcopy(m2)

    cruce1 = random.randint(0,len(m1.indices))


    for i in range(cruce1,len(m1.indices)):
        pivote = m1.indices[i]
        hijo1.indices[i] = m2.indices[i]
        hijo2.indices[i] = pivote

    return hijo1,hijo2


def roulette(popl):
    ganadores = []
    for p in popl:
        compara = random.random()
        if (p.fitness > compara):
            ganadores.append(p)
    return ganadores


datos = cargaSecuencias(archivo)
carga = int(sys.argv[2])
poplocal = []
print("Iniciando Trabajo")
#individuo = generaMatrizInicialNuevo(largo)
individuo = cargaIndividuo("results/PS00010.faarray6-72.fts")
AnnealingFitness(individuo)
print("inicio",individuo.fitness)
print(individuo.indices)
tinicial = 1
valor = 0
fitnessmejor = 0
mejorglobal = None
while(tinicial>0.00005):

    for i in range(10000):
        #print(tinicial,i,individuo.fitness)
        vecino = mutacion(individuo)
        AnnealingFitness(vecino)
        #print (vecino.fitness)
        resta = vecino.fitness - individuo.fitness
        #
        if (individuo.fitness > fitnessmejor):
            fitnessmejor = individuo.fitness
            mejorglobal = copy.deepcopy(individuo)
            escribeindividuo(mejorglobal,tinicial,"results/annealing"+archivo+"individuo"+str(largo)+"-"+str(len(datos))+".fts")

        if (resta>0):
            individuo = copy.deepcopy(vecino)

        else:
            r=random.random()
            valor = math.exp((decimal.Decimal(-resta)/decimal.Decimal(tinicial)))

            #print(r,valor)
            if(r>valor):
                #print("PEOR")
                individuo = copy.deepcopy(vecino)
    print(individuo.fitness)
    escribeFitness(individuo.fitness,largo,"results/annealing"+archivo+"mejor"+str(largo)+".fts")
    tinicial = tinicial*0.99
