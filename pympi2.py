import random
import copy
import decimal
import sys
import uuid
import hashlib
import time
import numpy
from mpi4py import MPI
comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

archivo = "PS00010.fa"
poblacion = 50000
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
    f=open(archivo,'a')
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

def newFitness(mmatriz):
    localr = []
    largototal = len(datos)
    alfabeto=["A","C","E","D","G","F","I","H","K","M","L","N","Q","P","S","R","T","W","V","Y","X"]
    for matriz in mmatriz:
        ind = 0
        palabras = []
        fitness=0
        for k in matriz.indices:

            palabras.append(datos[ind][k:k+largo])
            ind = ind + 1

        for i in range(largo):
                counts = {"A":0,"C":0,"D":0,"E":0,"F":0,"G":0,"H":0,"I":0,"K":0,"L":0,"M":0,"N":0,"P":0,"Q":0,"R":0,"S":0,"T":0,"V":0,"W":0,"Y":0,"X":0}

                for l in palabras:
                    counts[l[i]]+=1
                #print(counts)

                valor = counts.values()
                maximo = max(valor)
                if(maximo == largototal):
                    maximo = maximo + largo*1.1
                fitness = fitness + maximo
                ceros = list(counts.values()).count(0)

                fitness = fitness + ceros



        propor = decimal.Decimal(decimal.Decimal(fitness)/(decimal.Decimal((len(datos)+(len(datos)-1)+largo)*largo)))

        matriz.setFitness(propor)
def mutacion(m):

    gen = random.randint(0,len(m.indices)-1)
    nuevogen = random.randint(0,len(datos[gen])-largo)
    m.indices[gen] = nuevogen


def parallelFitness(mmatriz):
    localr = []
    alfabeto=["A","C","E","D","G","F","I","H","K","M","L","N","Q","P","S","R","T","W","V","Y","X"]
    for matriz in mmatriz:
        fitness = 0
        palabra = ""
        nkmers = len(matriz.kmers)
        lkmers = len(matriz.kmers[0].adn)
        for k in range(len(matriz.kmers[0].adn)):
            pal = ""
            counts = {"A":0,"C":0,"D":0,"E":0,"F":0,"G":0,"H":0,"I":0,"K":0,"L":0,"M":0,"N":0,"P":0,"Q":0,"R":0,"S":0,"T":0,"V":0,"W":0,"Y":0,"X":0}
            for i in matriz.kmers:
                counts[i.adn[k]]+=1
                pal = pal + i.adn

            ho = hashlib.md5(pal.encode())
            valor = counts.values()


            maximo = max(valor)

            if (maximo == nkmers):
                maximo = maximo + lkmers
            fitness = fitness + maximo

            ceros = list(counts.values()).count(0)

            fitness = fitness + ceros



        propor = decimal.Decimal(decimal.Decimal(fitness)/(decimal.Decimal((len(datos)+lkmers+20)*lkmers)))
        matriz.setFitness(propor)
        matriz.setCode(ho.hexdigest())


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

def seleccion(popl):
    ganadores = []
    distintos = {}
    largo = len(popl)
    while(len(popl)>2):
        participante1 = random.randint(0,len(popl)-2)
        participante2 = random.randint(0,len(popl)-2)

        if (popl[participante1].fitness >= popl[participante2].fitness):
            #if (popl[participante1].code not in distintos):
            ganadores.append(popl[participante1])
            #    distintos[popl[participante1].code]=1



        else:
            #if (popl[participante2].code  not in distintos):
            ganadores.append(popl[participante2])
            #distintos[popl[participante1].code]=1

        popl.pop(participante1)
        popl.pop(participante2)

    return ganadores




datos = cargaSecuencias(archivo)
carga = int(sys.argv[2])
if (rank==0):
    print("carga "+str(carga))
    popglobal = []


    print("Iniciando Trabajo")
    #cargaSecuencias(archivo)

    #for i in range(1,size):
    #    comm.send(datos,dest=i,tag=10)

    for c in range(10000):
        start_time = time.process_time()
        if (carga==1):
            print("Recuperando archivo")
            popglobal = cargaRespaldo("recover/PS00010.fa629recover.fts")
            carga=0
        if(len(popglobal)>0):
            pop1 = []
            pop2 = []
            pop3 = []
            pop4 = []
            for p in popglobal:
                destino = random.randint(1,4)
                if (destino==1):
                    pop1.append(p)
                if (destino==2):
                    pop2.append(p)
                if (destino==3):
                    pop3.append(p)
                if (destino==4):
                    pop4.append(p)

            comm.send(pop1,dest=1,tag=20)
            comm.send(pop2,dest=2,tag=20)
            comm.send(pop3,dest=3,tag=20)
            comm.send(pop4,dest=4,tag=20)
        else:
            for i in range(1,size):
                comm.send(popglobal,dest=i,tag=20)
        popglobal=[]
        for i in range (1,size):

            lista = comm.recv(source=i,tag=30)

            for l in lista:
                popglobal.append(l)
        suma = 0
        mejor = None
        mejorfitness = 0
        for p in popglobal:
            suma = suma + p.fitness
            if (p.fitness > mejorfitness):
                mejorfitness = p.fitness
                mejorg = copy.deepcopy(p)
        print(c,suma/len(popglobal),len(popglobal))
        if (c % 10 == 0):
            guardaGeneracion("recover/"+archivo+str(largo)+str(len(datos))+"recover.fts",popglobal)

        escribeFitness(suma/len(popglobal),largo,"results/"+archivo+str(largo)+"-"+str(len(datos))+"-fitnesspromedio.fts")
        #escribeindividuo(mejorg,c,"results/"+archivo+str(largo)+"-"+str(len(datos))+"-generacion.fts")
        escribeFitness(mejorg.fitness,largo,"results/"+archivo+str(largo)+"-"+str(len(datos))+"-mejorfitness.fts")

        end_time = time.process_time()
        print("Global Time")
        print(end_time-start_time)

if (rank!=0):

    while(True):


        poplocal = comm.recv(source=0,tag=20)

        if (len(poplocal)==0):
            poplocal = []
            for i in range(poblacion):
                m = generaMatrizInicialNuevo(largo)
                poplocal.append(m)

        for mu in range (int(len(poplocal)*0.1)):
            amutar = random.randint(0,len(poplocal)-1)
            mutacion(poplocal[amutar])

        newFitness(poplocal)
        #seleccionados = seleccion(poplocal)
        seleccionados = roulette(poplocal)
        print(len(seleccionados),rank)
        poplocal = []
        poplocal = seleccionados

        hijos = []
        for i in range(poblacion):
            cruza1 = random.randint(0,len(poplocal)-1)
            cruza2 = random.randint(0,len(poplocal)-1)

            hijo1,hijo2 = cruzamiento(poplocal[cruza1],poplocal[cruza2])
            hijos.append(hijo1)
            hijos.append(hijo2)
        poplocal = []
        poplocal = hijos

        comm.send(poplocal,dest=0,tag=30)
