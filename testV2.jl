## Benemérita Universidad Autónoma de Puebla, México.
### Julia V1.5.0
#### Dr. Enrique Ricardo Pablo Buendia Lozada

#Estadística para la validación de Pruebas en:
#- Cultura Física, 
#- Educación Física,
#- Fisioterapia, 
#- Recreación, 
#- Ciencias de la actividad Física, 
#- Deportes, 
#- entre otras.

##### Versión: V1 Julia 1.5.0
##### Caso: Paramétrico
##### Febrero 2021

using Gadfly, Cairo
using Statistics, Distributions
using CSV, DataFrames
using SimpleANOVA


## Bland - Altman

function bland_altman(x,y,i,j)
    
    n1=length(x)
    n2=length(y)
    if n1==n2 
        mn=0.5*(x+y)
        dd=x-y
        md=sum(dd)./length(dd)
       
        vd=var(dd)
        
        ls=md+2*(vd)^(0.5)
        li=md-2*(vd)^(0.5)
        
        # View Blant - Altman
        p=plot(x=mn,y=dd,Geom.point,yintercept=[ls, li],Geom.hline(color=["red","red"]),Guide.title("Bland-Altman"))
        img = SVG("img_muestras_$i-$j.svg", 6inch, 4inch)
        draw(img, p)
        draw(PNG("img_muestras_$i-$j.png", 6inch, 4inch), p)
        
    else
        println("Las muestras no son del mismo tamaño...")
        
    end
end

function ver_outliers(x,v1)    
    outl=0
    lsbox1=quantile(x,0.75)+1.5(quantile(x,.75)-quantile(x,.25))
    libox1=quantile(x,0.25)-1.5(quantile(x,.75)-quantile(x,.25))
    # a[a .> 4]
    os=x[x .> lsbox1]
    if length(os)>=1
        println("Outliers: ",os," en muestra $v1 ")
        outl=1
    end
    os=[]
    oi=x[x .< libox1]
    if length(oi)>=1
        println("Outliers: ",oi," en muestra $v1 ")
        outl=1
    end
    oi=[]
    if outl==0
        println("Sin outliers ...muestra $v1 ")
    end
end


# coeficiente de de correlación de concordancia de Lawrence - Lin
#
function ll(x,y)
    rc=(2*cov(x,y))/(var(x)+var(y)+(mean(x)-mean(y))^2)
    print("Lawrence - Lin $rc ")
    
    if rc >= 0. && rc < 0.01
        println("pobre --> Interpretación de Landis y Koch (1977)")
    elseif rc >=0.01 && rc<=0.20
        println("leve  --> Interpretación de Landis y Koch (1977)")
    elseif rc>=0.21 && rc <=0.40
        println("regular  --> Interpretación de Landis y Koch (1977)")
    elseif rc>=0.41 && rc<=0.60
        println("moderado  --> Interpretación de Landis y Koch (1977)")
    elseif rc>=0.61 && rc<=0.80
        println("substancial  --> Interpretación de Landis y Koch (1977)")
    elseif rc>=0.81 && rc<=1
        println("casi perfecto  --> Interpretación de Landis y Koch (1977)")
    end
end

#tabla de valoración
function tablaval(x)
    dest=(var(x))^(0.5)
    media=mean(x)
    dest1=media+dest
    dest2=media+dest*2
    dest3=media+dest*3
    ndest1=media-dest
    ndest2=media-dest*2
    ndest3=media-dest*3
     ff=x[x .> dest3] 
    f=length(ff)
    println(" $dest3 - o más     $f")
     ff=x[(x .> dest2) .& (x .<= dest3)] 
    f=length(ff)
    println(" $dest2 - $dest3     $f")
     ff=x[(x .> dest1) .& (x .<= dest2)] 
    f=length(ff)
    println(" $dest1 - $dest2     $f")
     ff=x[(x .> media) .& (x .<= dest1)] 
    f=length(ff)
    println(" $media - $dest1     $f")
    
     ff=x[(x .> ndest1) .& (x .<= media)]
    f=length(ff)
    println(" $ndest1 - $media     $f")
     ff=x[(x .> ndest2) .& (x .<= ndest1)] 
    f=length(ff)
    println(" $ndest2 - $ndest1     $f")
     ff=x[(x .> ndest3) .& (x .<= ndest2)]
    f=length(ff)
    println(" $ndest3 - $ndest2     $f")
     ff=x[x .<= ndest3]
    f=length(ff)
    println(" o menos - $ndest3     $f")
end

function normalk2(x,i)
# D'Agostino-Pearson's K2 test for assessing normality of data using skewness and kurtosis.
# Adaptado y traducido de http://welles.dm.unibo.it/~simoncin/DagosPtest.m
alpha = 0.05
n = length(x)
s1=sum(x)
s2 = sum(x.^2)
s3 = sum(x.^3)
s4 = sum(x.^4)
ss = s2-(s1^2/n)
v = ss/(n-1)
k3 = ((n*s3)-(3*s1*s2)+((2*(s1^3))/n))/((n-1)*(n-2))
g1 = k3/sqrt(v^3)
k4 = ((n+1)*((n*s4)-(4*s1*s3)+(6*(s1^2)*(s2/n))-((3*(s1^4))/(n^2)))/((n-1)*(n-2)*(n-3)))-((3*(ss^2))/((n-2)*(n-3)))
g2 = k4/v^2
eg1 = ((n-2)*g1)/sqrt(n*(n-1))  #measure of skewness
eg2 = ((n-2)*(n-3)*g2)/((n+1)*(n-1))+((3*(n-1))/(n+1))  #measure of kurtosis
A = eg1*sqrt(((n+1)*(n+3))/(6*(n-2)))
B = (3*((n^2)+(27*n)-70)*((n+1)*(n+3)))/((n-2)*(n+5)*(n+7)*(n+9))
C = sqrt(2*(B-1))-1
D = sqrt(C)
E = 1/sqrt(log(D))
F = A/sqrt(2/(C-1))
Zg1 = E*log(F+sqrt(F^2+1))

G = (24*n*(n-2)*(n-3))/((n+1)^2*(n+3)*(n+5))
H = ((n-2)*(n-3)*abs(g2))/((n+1)*(n-1)*sqrt(G))
J = ((6*(n^2-(5*n)+2))/((n+7)*(n+9)))*sqrt((6*(n+3)*(n+5))/((n*(n-2)*(n-3))))
K = 6+((8/J)*((2/J)+sqrt(1+(4/J^2))))
L = (1-(2/K))/(1+H*sqrt(2/(K-4)))
Zg2 = (1-(2/(9*K))-L^(1/3))/sqrt(2/(9*K))

K2 = Zg1^2 + Zg2^2  #D'Agostino-Pearson statistic
X2 = K2  #approximation to chi-distribution

df = 2.  #degrees of freedom
P=1-ccdf(Chisq(df), X2)
#println("P = $P ")
if P>alpha
    println("La muestra $i se distribuye normalmente ...")
    println("Tabla de Valoración    [NORMA]")
    println("Intervalos                            frecuencia")
    tablaval(x)
        println()
        println()
else
    println("La muestra $i NO viene de una distribución normal")
end
end








df=CSV.read("datos.csv";header=false)      # Ejemplo

#names(df)
#rename!(df,["c1","c2","c3"])              # se puede omitir, pues no importan los títulos

# realizar todas las combinaciones de 2
# de todas las muestras disponibles

muest   =ncol(df)
muestm1 = muest-1
i=1

while i<=muestm1
    j=i+1
    while j<=muest
       x = df[i]
       y = df[j]
       
        ver_outliers(x,i)
        ver_outliers(y,j)
        
	co=cor(x,y) # coeficiente de correlación de PEARSON
        print("Pearson: $co ")
        if co >= 0. && co < 0.01
            println("pobre --> Interpretación de Landis y Koch (1977)")
        elseif co >=0.01 && co<=0.20
            println("leve  --> Interpretación de Landis y Koch (1977)")
        elseif co >=0.21 && co <=0.40
            println("regular  --> Interpretación de Landis y Koch (1977)")
        elseif co >=0.41 && co <=0.60
            println("moderado  --> Interpretación de Landis y Koch (1977)")
        elseif co >=0.61 && co <=0.80
            println("substancial  --> Interpretación de Landis y Koch (1977)")
        elseif co >=0.81 && co <=1
            println("casi perfecto  --> Interpretación de Landis y Koch (1977)")
        end

        bland_altman(x,y,i,j) # guarda las gráficas para verlas en el navegador
        ll(x,y)   # coeficiente de correlación de  concordancia de Lawrence - Lin
        println(" imagen de muestras $i $j guardada ...")
        
        j=j+1
    end
    i=i+1
end 

# comprobar si todas las muestras provienen de la distribución normal
ccc=ncol(df)
i=1

nren=length(df[1])
while i<=ccc
     x = df[i]
     normalk2(x,i)
     i=i+1
end


mat1=convert(Matrix,df)
# verificar igualdad de varianzas [Homocedasticidad]
l=levene(mat1)
println("Si F > p rechazar Ho: Las muestras tienen varianzas iguales ...")
# consultar : https://github.com/BioTurboNick/SimpleANOVA.jl
println(l)

# verificar igualdad de promedios  [ANOVA de una vía]
p=anova(mat1)
println("Si F > p rechazar Ho: Las muestras tienen medias iguales ...")
# consultar : https://github.com/BioTurboNick/SimpleANOVA.jl
println(p)


