class AdnRecord:
    def __init__(self,kmer,position,largo,inicio):
        self.kmer = kmer
        self.position = position
        self.largo = largo
        self.inicio = inicio

        

        
def readGenome(filename):
    genome = ''
    with open(filename,'r') as f:
        for line in f:
            if not line[0] == '>':
                genome += line.rstrip()
                
    return genome        


def HammingDistance(kmer1,kmer2):
    contador = 0
    for i in range(len(kmer1)):
        if (kmer1[i]!=kmer2[i]):
            contador += 1
    return contador

def buscaVecino(palabra,adn,lmotif):
    
    mejor = 1000
    mejorpalabra = ""
    for i in range(len(adn)-lmotif):
        buscado = ""
        
        for l in range(i,i+lmotif):
            buscado = buscado+adn[l]
        distancia = HammingDistance(palabra,buscado)
        if (distancia<mejor):
            mejor=distancia
            mejorpalabra = buscado
    return mejorpalabra
    
def evaluaMatriz(datos):
    dato = [0] * 4
    palabratotal = ""
    alfabeto = ["A","C","G","T"]
    score = 0
    for i in range(len(datos[0])):
        counts = {"A":0,"C":0,"G":0,"T":0}
        for j in datos:
            
            counts[j[i]]+=1
            #print j[i]
            #print counts
            dato[0] = counts["A"]
            dato[1] = counts["C"]
            dato[2] = counts["G"]
            dato[3] = counts["T"]
        
        
            valor = max(dato)
            posicion = dato.index(valor)
        score = score + valor
        palabratotal = palabratotal + alfabeto[posicion]
#print palabratotal +" "+str(score)
    return palabratotal,score
        
        
            

def MotifSearch(datos,lmotif):
    motif = ""
    maxscore = 0
    palabras = [""] * (len(datos))
   
    for i in range(len(datos[0])-lmotif):
        palabra = ""
        for j in range(i,i+lmotif):
            
            palabra=palabra+datos[0][j] 
            ind = 1
        for k in datos[1:]:
            
            better = buscaVecino(palabra,k,lmotif)
            #print palabra+"-"+better
            palabras[ind] = better
            ind = ind + 1
        palabras[0] = palabra   
        #print palabras
        palabratotal,score = evaluaMatriz(palabras)
        #print maxscore
        if (score>maxscore):
            motif = palabratotal
            maxscore = score
        print palabratotal+" score "+str(score)
    print "MOTIF FINAL: "+motif+" "+str(maxscore)
        #for h in range(len(palabra)):
            #posicion,valor = evaluaMatriz(palabras)
          
            #palabratotal,score = evaluaMatriz(palabras)
            #print palabratotal+"-"+str(score)
           
            
       
        
    
    


def patternSearch():
    import mysql.connector
    cnx = mysql.connector.connect(user="root",password="barton",host="localhost",database="adn",port="3307")
    cursor = cnx.cursor()   
    cursor.execute("select max(position) from kmers")
    data = cursor.fetchall()
    for i in data:
        largo = i[0]
    print largo
    cursor.close()
    cont = 0
    for i in range(largo):
        
        cursor = cnx.cursor()
        cursor.execute("select * from kmers where position="+str(i))
        data = cursor.fetchall()
        for j in data:
            
            get = data[0]
            kmer = get[1]
        cursor.close()
        
        cursor = cnx.cursor()
        cursor.execute("select * from kmers where kmer='"+kmer+"'")
        data = cursor.fetchall()
        
        print str(cont)+"-"+kmer +"-"+ str(len(data))
        cont = cont + 1
            
        
    cnx.close()
    
    
        
        

        
        
def lookPattern():
    datos = []
    kmers = []
    import mysql.connector
    cnx = mysql.connector.connect(user="root",password="barton",host="localhost",database="adn",port="3307")
    cursor = cnx.cursor()
    cursor.execute("SELECT distinct(kmer),count(*) from kmers group by kmer")
    contador = 0
    data = cursor.fetchall()
    for i in data:
        get = i
        kmers.append(get[0])    
        datos.append(get[1])
        
        contador = contador + 1
    
    import numpy as np
    #import matplotlib.pyplot as plt
    
    print datos
    N = contador
    ind = np.arange(N)
    fig,ax = plt.subplots()
    rects = ax.bar(ind,datos,0.35,color='y')
    
    
    
    
    
    
    #plt.bar(range(len(datos)),datos)
    #plt.show()
        
        
def generaIndex(string,k):
    import mysql.connector
    cnx = mysql.connector.connect(user="root",password="barton",host="localhost",database="adn",port="3307")
    cursor = cnx.cursor()
    
    palabras = []
    finales = []
    for i in range(len(string)-k):
        palabra = ""
        for j in range(i,i+k):
            palabra = palabra+string[j]
        kmer = AdnRecord(palabra,i,len(palabra),palabra[0])
            
        add_adn = ("INSERT INTO kmers (kmer,position,largo,inicio) values(%(kmer)s,%(position)s,%(largo)s,%(inicio)s)")

        data_adn = {
            'kmer': kmer.kmer,
            'position': kmer.position,
            'largo': kmer.largo,
                'inicio': kmer.inicio
        }
        cursor.execute(add_adn,data_adn)
        cnx.commit()
    
            
    
   
  
    
    
    
       
        
        
genomep = readGenome("lambda_virus.fa")
genome = "ACTCTAGCTAGCTAGCTAGCATGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCATCGATCGATCGCGCTAGCATGCTAGCATGCTAGCTAGCATGCTAGCTAGCTAGCTAGCTACGATCGTACGATTCGATCGTCGTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCACTCTAGCTAGCTAGCTAGCATGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCATCGATCGATCGCGCTAGCATGCTAGCATGCTAGCTAGCATGCTAGCTAGCTAGCTAGCTACGATCGTACGATTCGATCGTCGTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCACTCTAGCTAGCTAGCTAGCATGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCATCGATCGATCGCGCTAGCATGCTAGCATGCTAGCTAGCATGCTAGCTAGCTAGCTAGCTACGATCGTACGATTCGATCGTCGTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCACTCTAGCTAGCTAGCTAGCATGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCATCGATCGATCGCGCTAGCATGCTAGCATGCTAGCTAGCATGCTAGCTAGCTAGCTAGCTACGATCGTACGATTCGATCGTCGTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCACTCTAGCTAGCTAGCTAGCATGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCATCGATCGATCGCGCTAGCATGCTAGCATGCTAGCTAGCATGCTAGCTAGCTAGCTAGCTACGATCGTACGATTCGATCGTCGTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCACTCTAGCTAGCTAGCTAGCATGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCATCGATCGATCGCGCTAGCATGCTAGCATGCTAGCTAGCATGCTAGCTAGCTAGCTAGCTACGATCGTACGATTCGATCGTCGTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCACTCTAGCTAGCTAGCTAGCATGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCATCGATCGATCGCGCTAGCATGCTAGCATGCTAGCTAGCATGCTAGCTAGCTAGCTAGCTACGATCGTACGATTCGATCGTCGTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCACTCTAGCTAGCTAGCTAGCATGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCATCGATCGATCGCGCTAGCATGCTAGCATGCTAGCTAGCATGCTAGCTAGCTAGCTAGCTACGATCGTACGATTCGATCGTCGTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGC"

#print len(genomep)
#generaIndex(genomep,11)
#lookPattern()
#HammingDistance("ATCG","ATAG")
datos = []
#adn1 = "cctgatagacgctatctggctatccacgtacataggtcctctgtgcgaatctatgcgtttccaaccat"
#adn2 = "agtactggtgtacatttgatacgtacgtacaccggcaacctgaaacaaacgctcagaaccagaagtgc"
#adn3 = "aaaagtccgtgcaccctctttcttcgtggctctggccaacgagggctgatgtataagacgaaaatttt"
#adn4 = "agcctccgatgtaagtcatagctgtaactattacctgccacccctattacatcttacgtacgtataca"
#adn5 = "ctgttatacaacgcgtcatggcggggtatgcgttttggtcgtcgtacgctcgatcgttaacgtaggtc"

#adn1 = "TCTCATCCGGTGGGAATCACTGCCGCATTTGGAGCATAAACAATGGGGGG"
#adn2 = "TACGAAGGACAAACACTTTAGAGGTAATGGAAACACAACCGGCGCATAAA"
#adn3 = "ATACAAACGAAAGCGAGAAGCTCGCAGAAGCATGGGAGTGTAAATAAGTG"
#adn4 = "GGCGCCTCATTCTCGGTTTATAAGCCAAAACCTTGTCGAGGCAACTGTCA"
#adn5 = "TCAAATGATGCTAGCCGTCGGAATCTGGCGAGTGCATAAAAAGAGTCAAC"


adn1 = "CGGATGGAATCGCCGCTTTTGAATTCACCTCCGGGGTATTATTATTATTCTTAGTAGTCGCGGTCGTGCGGACACCCGGAGTTATGCGGGCCCGAAAGCTCATTATGTAGTAAAGCTAGGTAATGTTAAGGGCGTAAGAGCCAACGCAAGGCAGCAATAGCCTGGTATTCCCACATATCAAGAAAGCTTAAAAAGTTGAGACAGGGAATTTGAAGGCGAAGATTGCCGAACTGGCCAATACCCACTACTTTTTTTTTGGTTTGCTTGGTTTCTTCCTGTCGCTTGCCAACTTGTGGCATCTTCCCCACACTATATTATAAGGATCGTCCTATGTATAGGCAATATTATCCATTTCACTCGCTAACAAATGTACGTATATATATGGAGCAACAAGTAGTGCAATTACAGACGTGTATTTTGTCTTGATCTTGCTTTTTGTATGATAGGCCTAAGAATAACAGTGCGAACATATAAGAAACATCCCTCATACTACCACACAT"
adn2 = "TCAATTAGTTCTGTTGCGCTTGACAATATATGTCGTGTAATACCGTCCCTTAGCAGAAGAAAGAAAGACGGATCCATATATGTTAAAATGCTTCAGAGATGTTTCTTTAATGTGCCGTCCAACAAAGGTATCTTCTGTAGCTTCCTCTATTTTCGATCAGATCTCATAGTGAGAAGGCGCAATTCAGTAGTTAAAAGCGGGGAACAGTGTGAATCCGGAGACGGCAAGATTGCCCGGCCCTTTTTGCGGAAAAGATAAAACAAGATATATTGCACTTTTTCCACCAAGAAAAACAGGAAGTGGATTAAAAAATCAACAAAGTATAACGCCTATTGTCCCAATAAGCGTCGGTTGTTCTTCTTTATTATTTTACCAAGTACGCTCGAGGGTACATTCTAATGCATTAAAAGAC"
adn3 = "TGGAACTTTCAGTAATACGCTTAACTGCTCATTGCTATATTGAAGTACGGATTAGAAGCCGCCGAGCGGGCGACAGCCCTCCGACGGAAGACTCTCCTCCGTGCGTCCTCGTCTTCACCGGTCGCGTTCCTGAAACGCAGATGTGCCTCGCGCCGCACTGCTCCGAACAATAAAGATTCTACAATACTAGCTTTTATGGTTATGAAGAGGAAAAATTGGCAGTAACCTGGCCCCACAAACCTTCAAATTAACGAATCAAATTAACAACCATAGGATGATAATGCGATTAGTTTTTTAGCCTTATTTCTGGGGTAATTAATCAGCGAAGCGATGATTTTTGATCTATTAACAGATATATAAATGGAAAAGCTGCATAACCACTTTAACTAATACTTTCAACATTTTCAGTTTGTATTACTTCTTATTCAAATGTCATAAAAGTATCAACAAAAAATTGTTAATATACCTCTATACTTTAACGTCAAGGAGAAAAAACTATA"
adn4 = "GCTTTTACTATTATCTTCTACGCTGACAGTAATATCAAACAGTGACACATATTAAACACAGTGGTTTCTTTGCATAAACACCATCAGCCTCAAGTCGTCAAGTAAAGATTTCGTGTTCATGCAGATAGATAACAATCTATATGTTGATAATTAGCGTTGCCTCATCAATGCGAGATCCGTTTAACCGGACCCTAGTGCACTTACCCCACGTTCGGTCCACTGTGTGCCGAACATGCTCCTTCACTATTTTAACATGTGGAATTCTTGAAAGAATGAAATCGCCATGCCAAGCCATCACACGGTCTTTTATGCAATTGATTGACCGCCTGCAACACATAGGCAGTAAAATTTTTACTGAAACGTATATAATCATCATAAGCGACAAGTGAGGCAACACCTTTGTTACCACATTGACAACCCCAGGTATTCATACTTCCTATTAGCGGAATCAGGAGTGCAAAAAGAGAAAATAAAAGTAAAAAGGTAGGGCAACACATAGT"
adn5 = "GGACCCTGACGGCGACACAGAGATGACAGACGGTGGCGCAGGATCCGGTTTAAACGAGGATCCCTTAAGTTTAAACAACAACAGCAAGCAGGTGTGCAAGACACTAGAGACTCCTAACATGATGTATGCCAATAAAACACAAGAGATAAACAACATTGCATGGAGGCCCCAGAGGGGCGATTGGTTTGGGTGCGTGAGCGGCAAGAAGTTTCAAAACGTCCGCGTCCTTTGAGACAGCATTCGCCCAGTATTTTTTTTATTCTACAAACCTTCTATAATTTCAAAGTATTTACATAATTCTGTATCAGTTTAATCACCATAATATCGTTTTCTTTGTTTAGTGCAATTAATTTTTCCTATTGTTACTTCGGGCCTTTTTCTGTTTTATGAGCTATTTTTTCCGTCATCCTTCCCCAGATTTTCAGCTTCATCTCCAGATTGTGTCTACGTAATGCACGCCATCATTTTAAGAGAGGACAGAGAAGCAAGCCTCCTGAAAG"
adn6 = "ATCGATTTTGCAGATTGTTCTAAAAGTAAATGGATTGCTATTTTCTTTCCGAGACTACTCTAAAAAAATTTATTGAGTATGAGATCGTTTTTAGATAAATTATATATATTGTAAAGCTATTAACTAATCTCCTATATCAATTTCTTCTTGCTTAACCCCGTGTGGTTGTTTAGGTCCATCTCCTTTTTCCTTTTAATTTTTTTACCTTTATTAATTCCTTCACCTCTCTAAACCCCAGTTTTATATCGTATATGCTATCTACAGGTCCACTTTACACTTAATAATATAAAAATACTACTATAAAGGAACCAGAAAAATAAAAAAGGGTCATTATTTATTTGAGCAGATCATTATCAAACGCATAGGAAGAGAAAAAACACAGTTTTATTTTTTTTCCACACATATTTATTGGTCTCCTAGTACATCAAAGAGCATTTTAATGGGTTGCTGATTTGTTTTACCTACATTTTCTAGTACAAAAAAAAAACAAAAAAAGAATC"
adn7 = "GCAGCTTCACTTTTAAGTTTCTTTTTCTCCTCACGGCGCAACCGCTAACTTAAGCTAATCCTTATGAATCCGGAGAAAAGCGGGGTCTTTTAACTCAATAAAATTTTCCGAAATCCTTTTTCCTACGCGTTTTCTTCGGGAACTAGATAGGTGGCTCTTCCACCTGTTTTTCCATCATTTTAGTTTTTCGCAAGCCATGCGTGCCTTTTCGTTTTTGCGATGGCGAAGCAGGGCTGGAAAAATTAACGGTACGCCGCCTAACGATAGTAATAGGCCACGCAACTGGCGTGGACGACAACAATAAGTCGCCCATTTTTTATGTTTTCAAAACCTAGCAACCCCCACCAAACTTGTCATCGTTCCCGGATTCACAAATGATATAAAAAGCGATTACAATTCTACATTCTAACCAGATTTGAGATTTCCTCTTTCTCAATTCCTCTTATATTAGATTATAAGAACAACAAATTAAATTACAAAAAGACTTATAAAGCAACATA"
adn8 = "CCCCCCGGTTTTTTTTCCATGGGGCCCCATATTCCCCCGCCTGCAGGAAAACTTGGGGAAAGAGGAAAAACACTTCGGATAAAAACGGTCAAGAAGCTCTTCGACGATTTAGTGCCACCTTCATGAAAAATTCCAGAGTTTTTTCCAGCTGCTTTGATTTTACAGTCCATTATTCGGCGTCTAACGATTCTGATTAAGAAACAACGGAGGAAAACTCAAATTCTAATATAATATTTTTAAGTTTATGAAGGTGGGGTGGTAAGAAAAGCAACTAAAATAATCTTCAAGTCAATTAGTGGTGAAAAGCTTCAACACTGGGGAATGAATAATATGTCATCTAGAAAAAATTTTATATAAATACTCAGTGTTTTATTCATTATTCTCGATTCATTCACTTCAATTCCTCTTCATGAGTAATAGAAACCATAAAGAAAAGATATATTCAAAGCCTCTTATCAAGGTTTGGTTTTGAAACACTTTTACAATAAAATCTGCCAAAA"
adn9 = "TAAGGAGAAATATTATCTAAAAGCGAGAGTTTAAGCGAGTTGCAAGAATCTCTACGGTACAGATGCAACTTACTATAGCCAAGGTCTATTCGTATTGGTATCCAAGCAGTGAAGCTACTCAGGGGAAAACATATTTTCAGAGATCAAAGTTATGTCAGTTTCTTTTTCATGTGTAACTTAACGTTTGTGCAGGTATCATACCGGCCTCCACATAATTTTTGTGGGGAAGACGTTGTTGTAGCAGTCTCCTTATACTCTCCAACAGGTGCCTAAAGACTTCTTCAGGCCTCATAGTCTACATCTGGAGACAACATTAGATAGAAGTTTCCACAGAGGCAGCTTTCAATATACTTTCGGCTGTGTACATTTCATCCTGAGTGAGCGCATATTGCATAAGTACTCAGCATATAAAGAGACACAATATACTCCATACTTGTTGTGAGTGGTTTTAGCGTATTCAGTATAACAATAAGAATTACATCCAAGACTATTAATTAACT"









datos.append(adn8)
datos.append(adn9)
datos.append(adn1)
datos.append(adn2)
datos.append(adn3)
datos.append(adn4)
datos.append(adn5)
datos.append(adn6)
datos.append(adn7)
#readSeqs("MIG1-seqs")
MotifSearch(datos,10)