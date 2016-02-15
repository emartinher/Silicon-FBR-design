# -*- coding: utf-8 -*-
"""
Created on Wed Dec 30 20:01:39 2015

@author: emartinher

Diseño del reactor de lecho fluidizado
"""

#Importar librerias
import numpy as np
from scipy import optimize
from scipy import integrate
import matplotlib.pyplot as plt



# PARAMETROS (Sistema Internacional de unidades)

# Parametros generales
g=9.81
R=8.314 
T=973 #Kelvin
P=136325 #Pa (Parker, 2011), presion del reactor + presion atmosferica

# Flujo volumetrico de los gases
#Flujo molar de SiH4=13.9210e3 mol/h
SiH4_n38=13.9210e3 #mol/h
SiH4_v38=((SiH4_n38*R*T)/P)/3600 #m3/s
ypsilon_0=(SiH4_v38/0.2)
H2_v39=ypsilon_0*0.8
H2_n39=((P*H2_v39)/(R*T))*3600 #mol/h
print('H2_n39=', H2_n39)

# Caracterización del distribuidor
N_or=2575  #(orificios/m2) Numero de agujeros del distribuidor

# Caracterizacion de las particulas 
d_p=0.00025 #250 micrometros #diametro de la esfera cuyo volumen es igual al de la particula
fi_s=0.6 #esfericidad (Fogler 2005)
#d_eff= #diametro efectivo
rho_s=2330 #densidad de partícula

# Caracterización del gas reactante (SiH4)
Pm_g_r= 32.117 #Peso molecular del gas reactante
rho_g_r=((P*Pm_g_r)/(R*T))/1000 #densidad del gas reactante. Se considera gas ideal (P=5 psig (Parker, 2011), T=1073K)  #densidad del gas reactante
y_g_r=0.2 #fraccion molar de gas reactante en la mezcla de gases (%)
C_s_r=155 #Coeficiente de Ecuacion de Sutherland del gas reactante (Minkina, 1995)
mi_0_g_r=0.000010704 #Ns/m2. Viscosidad del gas reactante a 273K
mi_g_r=mi_0_g_r*((273+C_s_r)/(T+C_s_r))*((T/273)**1.5) #viscosidad del gas reactante. Ecuacion de Sutherland
#D_ab= #difusividad del silano en hidrogeno
print('rho_g_r (kg/m3)', rho_g_r)
print('mi_g_r (kg/(m·s)', mi_g_r)

# Caracterización del gas fluidizante (H2)
Pm_g_f= 2.02 #Peso molecular del gas fluidizante
rho_g_f=((P*Pm_g_f)/(R*T))/1000 #densidad del gas fluidizante. Se considera gas ideal (P=5 psig (Parker, 2011), T=1073K) 
y_g_f=0.8 #fraccion molar de gas reactante en la mezcla de gases (%)
C_s_f=84.4 #Coeficiente de Ecuacion de Sutherland del gas fluidizante (Minkina, 1995)
mi_0_g_f=8.3969e-6 #Ns/m2. Viscosidad del gas fluidizante a 273K
mi_g_f= mi_0_g_f*((273+C_s_f)/(T+C_s_f))*((T/273)**1.5) #viscosidad del gas fluidizante. Ecuacion de Sutherland 
print('rho_g_f (kg/m3)', rho_g_f)
print('mi_g_f (kg/(m·s)', mi_g_f)

#Caracterización de la mezcla de gases
rho_g= (rho_g_r*y_g_r)+(rho_g_f*y_g_f)
fi_r_f=((1+((mi_g_r/mi_g_f)**0.5)*((Pm_g_f/Pm_g_r)**0.25))**2)/((8*(1+(Pm_g_r/Pm_g_f)))**0.5)
fi_f_r=((1+((mi_g_f/mi_g_r)**0.5)*((Pm_g_r/Pm_g_f)**0.25))**2)/((8*(1+(Pm_g_f/Pm_g_r)))**0.5)
mi_g=((y_g_r*mi_g_r)/(y_g_r+y_g_f*fi_r_f))+((y_g_f*mi_g_f)/(y_g_f+y_g_r*fi_f_r)) #Viscosidad de la mezcla de gases. Ecuacion de Wilk
print('rho_g (kg/m3)', rho_g)
print('mi_g (kg/(m·s)', mi_g)
        


# MINIMA VELOCIDAD DE FLUIDIZACIÓN

#Ecuacion de minima velocidad de fluidización
#(1.75/(epsilon_mf**3*fi_s))*(((d_p*u_mf*rho_g)/mi_g)**2)+((150*(1-epsilon_mf))/(epsilon_mf**3*fi_s**2))*((d_p*u_mf*rho_g)/mi_g)=((d_p**3*rho_g*(rho_s-rho_g)*g)/mi_g**2)

#Calculo epsilon_mf
eta=g*(rho_s-rho_g)
epsilon_mf=(0.586*fi_s**(-0.72))*(((mi_g**2)/(rho_g*eta*d_p**3))**0.029)*((rho_g/rho_s)**0.021)
print('epsilon_mf=',epsilon_mf)

#Resolución de la ecuación para hallar la velocidad minima de fluidizacion
def  f_u_mf(x_u_mf, epsilon_mf, fi_s, d_p, rho_g, mi_g):
    """ Resolución de la ecuación para hallar la velocidad minima de fluidizacion """
    return (1.75/(epsilon_mf**3*fi_s))*(((d_p*x_u_mf*rho_g)/mi_g)**2)+((150*(1-epsilon_mf))/(epsilon_mf**3*fi_s**2))*((d_p*x_u_mf*rho_g)/mi_g)-((d_p**3*rho_g*(rho_s-rho_g)*g)/mi_g**2)
    
u_mf=optimize.root(f_u_mf, 1, args=(epsilon_mf, fi_s, d_p, rho_g, mi_g))
u_mf=np.float_(u_mf.x)
print ('u_mf (m/s) =',u_mf)




# VELOCIDAD TERMINAL DE FLUIDIZACION

d_p_prima=d_p*((rho_g*(rho_s-rho_g)*g)/(mi_g**2))**(1/3)

u_t_prima=((18/d_p_prima**2)+((2.335-1.744*fi_s)/(d_p_prima**0.5)))**(-1)

# Ecuación de velocidad terminal
#u_t_prima=u_t*((rho_g**2)/(mi_g*(rho_s-rho_g)*g))**(1/3)

#Resolución de la ecuación para hallar la velocidad minima de fluidizacion
def  f_u_t(x_u_t, rho_g, mi_g, rho_s, u_t_prima):
    """ Resolución de la ecuación para hallar la velocidad minima de fluidizacion """
    return x_u_t*((rho_g**2)/(mi_g*(rho_s-rho_g)*g))**(1/3)-u_t_prima
    
u_t=optimize.root(f_u_t, 1, args=(rho_g, mi_g, rho_s, u_t_prima))
u_t=np.float_(u_t.x)
print ('u_t (m/s) =',u_t)
    



# VELOCIDAD SUPERFICIAL DEL GAS

# u_mf<u_0<u_t

u_0=5*u_mf #Parker 2011




# DIAMETRO DEL LECHO

#Area del lecho
A_c=ypsilon_0/u_0

#Diametro del lecho
d_t=np.sqrt((4*A_c)/(np.pi))




# DIAMETRO DE BURBUJA

#Diametro de burbuja maximo(esta ecuacion debe emplearse usando cm. Tras el calculo se transforma a m)
d_bm=(0.65*(((np.pi/4)*((d_t*100)**2)*((u_0*100)-(u_mf*100)))**0.4))/100
print ('d_bm (m) =', d_bm)

#Diametro de burbuja inicial(esta ecuacion debe emplearse usando cm. Tras el calculo se transforma a m)
d_b0=((1.3/((g*100)**0.2))*((((u_0*100)-(u_mf*100))/(N_or/(100**2)))**0.4))/100
print ('d_b0 (m) =', d_b0)





# BUCLE ITERATIVO

z=1 #☻ Estimacion incial para que entre en el bucle
d_b=d_bm-(np.exp((-0.3*z)/d_t))*(d_bm-d_b0)
print ('d_b (estimacion inicial) (m) =', d_b)
    
d_b_calculado=0
contador=0
L_f=None

while np.absolute(d_b-d_b_calculado)>0.001:
  
    #Diametro de burbuja
    if contador==0:    
        z=1  #Altura  la que se evalua el diametro de burbuja. Se evalua a la mitad de la altura del lecho Es un parametro que se tantea, debiendo darle un valor inicial
    else:
        z=(L_f/2)
    
    
    contador=contador+1    
    
    d_b=d_bm-(np.exp((-0.3*z)/d_t))*(d_bm-d_b0)
    print ('d_b (m) =', d_b)




    # VELOCIDAD DE ASCENSO DE BURBUJA

    #Velocidad de ascenso de burbuja respecto a la emulsion
    u_br=0.711*((g*d_b)**(1/2))

    #Velocidad de ascenso de burbuja a través del lecho. Esta ecuacion es valida para solidos tipos Geldart B
    u_b=1.6*((u_0-u_mf)+1.13*(d_b**0.5))*(d_t**1.35)+u_br
    #u_b=u_0-u_mf+u_br
    print ('u_b (m/s) =', u_b)

    #Velocidad total de ascenso de la burbuja
    u_b_total=u_b+u_mf
    print ('u_b_total (m/s) =', u_b_total)




    # FRACCION DE BURBUJAS EN EL LECHO
    delta=(u_0-u_mf)/u_b




    # COEFICIENTE DE INTERCAMBIO BURBUJA-EMULSION
    K_be=4.5*(u_mf/d_b)




    # FRACCION DE SOLIDOS EN LA BURBUJA
    gamma_b=0.005 #(Kimura y Kojima, 1991)




    # REACCION
    #Ecuacion (Parametros->(Kimura y Kojima, 1991))
    #dX/dt=(ks_0/(1+K_H*P_H+K_S*P_S))*(S/V)*(1-X)+(kv_0/(1+K_V*P_H))*(1-X)
    #P_H=(y_g_f*P)/1000 #kPa para introducirla en la ecuación
    #P_S=(y_g_r*P)/1000 #kPa para introducirla en la ecuación
    #ks_0=2.15e8*np.exp(-1.915e5/(R*T))
    #K_H=0.034
    #K_S=7.6e-3*np.exp(32900/(R*T))
    #kv_0=2.14e13*np.exp(-221300/(R*T))
    #K_V=0.50

    #RKS=ks_0/(1+K_H*P_H+K_S*P_S)
    #RKV=kv_0/(1+K_V*P_H)

    #Constante cinetica para la fase burbuja
    #K_r_bubble=((1-gamma_b)/gamma_b)*(RKS*(gamma_b/(1-gamma_b))*(6/d_p)+RKV)

    #Constante cinetica para la fase emulsion
    #K_r_emulsion=(epsilon_mf/(1-epsilon_mf))*(RKS*((1-epsilon_mf)/epsilon_mf)*(6/d_p)+RKV)

    #Ecuacion (Parametros->(Parker, 2011))
    #-dCs/dt=(k_het*(S_Si/V_R)+k_hom)*Cs

    k_het=2.793e6*np.exp(-1.954e4/T) #T (K), (m3 reactor/m2 Si surfaces s)
    k_hom=2e13*np.exp(-2.604e4/T) #T (K), (s-1)

    #BALANCES 
    #Balance a la fase burbuja (1)
    #-u_b_total*delta*dCs_b/dz=K_r_bubble*Cs_bubble*delta*gamma_b+K_be*(Cs_bubble-Cs_emulsion)*delta

    #Balance a la fase emulsion (2)
    #-(1-delta)*u_mf*dCs_e/dz=-K_be*(Cs_bubble-Cs_emulsion)*delta+K_r_emulsion*Cs_emulsion*(1-delta)*(1-epsilon_mf)

    #Reordenando ec 2 (3)
    #Cs_emulsion=(K_be/(K_r_emulsion*(Vs/Vg)+K_be))*Cs_bubble

    #Balances (Modelo Kunii-Levenspiel para Geldart B)
    #Balance a la fase burbuja
    #-u_b_total*delta*dCs_b/dz=(k_het*(S_Si/V_R)+k_hom)*Cs_bubble*delta*gamma_b+K_be*(Cs_bubble-Cs_emulsion)*delta 
    print('Balance a la fase burbuja')
    print('-u_b_total*delta*dCs_b/dz=(k_het*(S_Si/V_R)+k_hom)*Cs_bubble*delta*gamma_b+K_be*(Cs_bubble-Cs_emulsion)*delta')

    #Balance a la fase emulsion 
    #-(1-delta)*u_mf*dCs_e/dz=-K_be*(Cs_bubble-Cs_emulsion)*delta+(k_het*(S_Si/V_R)+k_hom)*Cs_emulsion*(1-delta)*(1-epsilon_mf)
    print ('Balance a la fase emulsion')
    print ('-(1-delta)*u_mf*dCs_e/dz=-K_be*(Cs_bubble-Cs_emulsion)*delta+(k_het*(S_Si/V_R)+k_hom)*Cs_emulsion*(1-delta)*(1-epsilon_mf)')

    print('Condiciones de contorno para integrar los balances: Cs_bubble=Cs_emulsion=Cs_entrada de gases')

    #Conversion global (podemos determinar la altura del lecho en condiciones de fluidizacion)
    print('Revisar!!!!')
    Cs_i=(SiH4_n38/3600)/(ypsilon_0) #mol/m3. Concentracion de SiH4 a la entrada del reactor    
    # 1-X_s=(delta*u_b_total*Cs_bubble0+(1-delta)*u_mf*Cs_emulsion0)/(u_0*Cs_i)

    print('Definimos una conversion total del silano que queremos alcanzar')
    X_s=0.8
    print('X_s=', X_s)
    
    #S_Si/V_R=S_V

    #S_V=(6/(d_p))*((1-epsilon_mf)/epsilon_mf) #d_p Caussat et al, 1995
    #S_V=((4*np.pi*(d_p/2)**2)*(1-epsilon_mf)*(z*2*np.pi*(d_t/2)**2))/(z*2*np.pi*(d_t/2)**2)
    #S_V=((1/((4/3)*np.pi*(d_p/2)**3))*(4*np.pi*(d_p/2)**2))*(1-epsilon_mf)
    S_V=(6/(d_p))*(1-epsilon_mf)

    fi=((((k_het*(S_V)+k_hom)/K_be)**2)*(((1-epsilon_mf)-gamma_b*(u_mf/u_b_total))**2)+(((delta/(1-delta))+(u_mf/u_b_total))**2)+2*((k_het*(S_V)+k_hom)/K_be)*((1-epsilon_mf)-gamma_b*(u_mf/u_b_total))*((delta/(1-delta))-(u_mf/u_b_total)))**0.5

    q_1=(1/2)*((k_het*(S_V)+k_hom)/u_mf)*((1-epsilon_mf)+gamma_b*(u_mf/u_b_total))+(1/2)*(K_be/u_mf)*(((delta/(1-delta))+(u_mf/u_b_total))-fi)
    q_2=(1/2)*((k_het*(S_V)+k_hom)/u_mf)*((1-epsilon_mf)+gamma_b*(u_mf/u_b_total))+(1/2)*(K_be/u_mf)*(((delta/(1-delta))+(u_mf/u_b_total))+fi)

    psi_1=(1/2)-(1/2)*((1-delta)/delta)*((u_mf/u_b_total)-((k_het*(S_V)+k_hom)/K_be)*((1-epsilon_mf)-gamma_b*(u_mf/u_b_total))-fi)
    psi_2=(1/2)-(1/2)*((1-delta)/delta)*((u_mf/u_b_total)-((k_het*(S_V)+k_hom)/K_be)*((1-epsilon_mf)-gamma_b*(u_mf/u_b_total))+fi)

    #Calcular L_f
    # 1-X_s=(delta/(1-delta))*(1/(u_0*fi))*((1-psi_2)*(psi_1*delta*u_b_total+(1-delta)*u_mf)*np.exp(-q_1*L_f)+(psi_1-1)*(psi_2*delta*u_b_total+(1-delta)*u_mf)*np.exp(-q_2*L_f))

    def  f_L_f(x_L_f, delta, u_0, fi, psi_2, psi_1, u_b_total, u_mf, q_1, q_2, X_s):
        """ Resolución de la ecuación para hallar la altura del lecho en condiciones de fluidizacion """
        return (delta/(1-delta))*(1/(u_0*fi))*((1-psi_2)*(psi_1*delta*u_b_total+(1-delta)*u_mf)*np.exp(-q_1*x_L_f)+(psi_1-1)*(psi_2*delta*u_b_total+(1-delta)*u_mf)*np.exp(-q_2*x_L_f))-(1-X_s)
        
    L_f=optimize.root(f_L_f, 1, args=(delta, u_0, fi, psi_2, psi_1, u_b_total, u_mf, q_1, q_2, X_s))
    L_f=np.float_(L_f.x)    
    print('Altura del lecho en condiciones de fluidizacion:')
    print ('L_f (m/s) =',L_f)




    #COMPARACION ENTRE EL RESULTADO OBTENIDO Y EL ENCONTRADO EN LA ITERACIÓN ANTERIOR
    d_b_calculado=d_bm-(np.exp((-0.3*(L_f/2))/d_t))*(d_bm-d_b0)
    print('Comparacion entre el resultado obtenido y el encontrado en la iteracion anterior:')
    print ('d_b_calculado (m) =', d_b_calculado)
    print ('d_b (m) =', d_b)

print ('Numero de iteraciones=',contador)

print('Diametro de burbuja definitivo:')
print ('d_b (m) =', d_b)

print('Altura del lecho en condiciones de fluidizacion definitivo:')
print ('L_f (m/s) =',L_f)




#PESO DEL LECHO DE PARTICULAS
W=rho_s*A_c*L_f*(1-epsilon_mf)*(1-delta)
print ('Peso de las particulas del lecho')
print ('W (kg) = ', W)




#TIEMPO DE RESIDENCIA
t=L_f/u_b_total
print ('Tiempo de residencia')
print ('t (s) = ', t)




#PERDIDA DE CARGA EN EL LECHO
#Ecuación de Ergun
delta_P_lecho=L_f*(150*(((1-epsilon_mf)**2)/(epsilon_mf**3))*((mi_g*u_0)/((fi_s*d_p)**2))+1.75*((1-epsilon_mf)/(epsilon_mf**3))*((rho_g*u_0**2)/(fi_s*d_p)))
print ('Perdida de carga a traves del lecho')
print ('delta_P_lecho (Pa) = ', delta_P_lecho)




#PERDIDA DE CARGA EN EL DISTRIBUIDOR
#Ecuación de Ergun
delta_P_distribuidor=0.3*delta_P_lecho
print ('Perdida de carga a traves del distribuidor')
print ('delta_P_distribuidor (Pa) = ', delta_P_distribuidor)




#PERFILES DE CONCENTRACION DE SILANO
#Ecuaciones diferenciales
def ec_cinetica_depos (Cs,t):
    """Cinética de deposición de silano mediante las rutas heterogenea y homogenea
    Inputs:
        Cs[0]: Composicion de silano inicial en las burbujas
        Cs[1]: Composicion de silano inial en la emulsion
        t: vector de tiempos
    Outputs:
        dCsdt[0]: Composicion en las burbujas
        dCsdt[1]: Composicion en la emulsion
        dC_sdt: Composición de silano con el tiempo
    Condiciones de contorno:
        Cs[0]=Cs[1]=Cs_i
        """
    dCsdt=np.zeros_like(Cs)
    dCsdt[0] = (delta*gamma_b*(k_hom)*Cs[0]+delta*K_be*(Cs[0]-Cs[1]))/(-delta*u_b_total)
    dCsdt[1] = ((1-delta)*(1-epsilon_mf)*(k_het*S_V)*Cs[1]-delta*K_be*(Cs[0]-Cs[1]))/(-(1-delta)*u_mf)
    #dCsdt=-(k_het*S_V+k_hom)*K_be*Cs
    
    return dCsdt
    
t=np.linspace(0,L_f,500)
#Cs_i=3.4520086600914124
Cs_ini=[Cs_i,Cs_i]

print('t=altura')
res_Cs=integrate.odeint(ec_cinetica_depos,Cs_ini,t)
#print ('res_Cs',res_Cs)

#fig,axes=plt.subplots(1,2,sharex=False,sharey=False)
#axes[0].plot(t,res_Cs[:,0],'b')
#axes[1].plot(t,res_Cs[:,1],'r')

fig,axes=plt.subplots()
axes.plot(t,res_Cs[:,0],'b')
axes.set_xlabel('Altura (m)')
axes.set_ylabel('Concentración SiH4 (mol/m3)')
axes.grid('on')
fig.savefig('Perfil_concentracion_burbuja.png')

fig,axes=plt.subplots()
axes.set_xlabel('Altura (m)')
axes.set_ylabel('Concentración SiH4 (mol/m3)')
axes.plot(t,res_Cs[:,1],'r',linewidth=3)
axes.grid('on')
fig.savefig('Perfil_concentracion_emulsion.png')



#SCAVENGING
#Formacion de finos

#def ec_cinetica_depos (Cs,h):
#    """Cinética de deposición de silano mediante las rutas heterogenea y homogenea
#    Inputs:
#        Cs[0]: Composicion de silano inicial en las burbujas
#        Cs[1]: Composicion de silano inial en la emulsion
#        t: vector de tiempos
#    Outputs:
#        dCsdt[0]: Composicion en las burbujas
#        dCsdt[1]: Composicion en la emulsion
#        dC_sdt: Composición de silano con el tiempo
#    Condiciones de contorno:
#        Cs[0]=Cs[1]=Cs_i
#        """
#    dCsidh=(k_hom*Cs_b)/u_b
#    
#    return dCsdt

#Csi_finos=(k_hom*res_Cs[:,0])/u_b_total #Variacion de la concentracion de finos de silicio con la altura
Csi_finos=Cs_i-res_Cs[:,0]




#CONCENTRACION DE FINOS A LO LARGO DEL REACTOR
alpha=2e-4
r_scav_finos=(alpha*S_V*(Csi_finos*28.09)/1000)/u_b_total #kg/(m3 m)

fig,axes=plt.subplots()
axes.set_xlabel('Altura (m)')
axes.set_ylabel('Concentración Si (kg/m3)')
axes.plot(t,r_scav_finos,'k',linewidth=1)
axes.grid('on')
fig.savefig('Perfil_concentracion_finos.png')
#Fraccion de silicio perdido en  los finos arrastrados; flujo molar gases de salida (balances de materia):58470 mol/h
finos_perdidos=(r_scav_finos[-1]*((58470*R*T)/P))/((SiH4_n38/1000)*28.09*X_s)
print('Fraccion de silicio perdido en  los finos arrastrados',finos_perdidos)



#VELOCIDAD DE CAMBIO DE MASA DE LAS PARTÍCULAS
print('Revisar la concentracion y rho_finos!!')
#r_het=-k_het*S_V*(((Cs_i+(Cs_i*0.2))/2)/1000) #kmol/(m3 s)
r_het_calc=(k_het*S_V*(res_Cs[:,1]/1000)) #kmol/(m3 s)
tiempo=t/u_mf

fig,axes=plt.subplots()
axes.plot(tiempo,r_het_calc,'r',linewidth=5)
axes.set_xlabel('Tiempo (s)')
axes.set_ylabel('Concentración Si (kmol/m3)')
axes.grid('on')
fig.savefig('Perfil_concentracion_deposhet.png')

r_het_media=np.sum(r_het_calc)/(500)



rho_finos=2330
#r_scav=-alpha*S_V*(rho_finos/28.09)/1000 #kmol/(m3 s)
r_scav=(alpha*S_V*(Csi_finos)/1000) #kmol/(m3 s)

fig,axes=plt.subplots()
axes.set_xlabel('Tiempo (s)')
axes.set_ylabel('Concentración Si (kmol/m3)')
axes.plot(tiempo,r_scav,'b',linewidth=1)
axes.grid('on')
fig.savefig('Perfil_concentracion_deposhom.png')

r_scav_media=np.sum(r_scav)/(500)


r_part=(r_het_calc+r_scav)*28.09*(1-epsilon_mf)*(1-delta)

r_part_media=(r_het_media+r_scav_media)*28.09*(1-epsilon_mf)*(1-delta) #kg/(m3 s)

fig,axes=plt.subplots()
axes.set_xlabel('Tiempo (s)')
axes.set_ylabel('Concentración Si (kg/m3)')
axes.plot(tiempo,r_part,'k',linewidth=1)
axes.grid('on')
fig.savefig('Perfil_concentracion_depostot.png')


#TIEMPO DE CRECIMIENTO DE UNA PARTICULA
#Fuerza de arrastre 
F_a=delta_P_lecho*(d_t/2)**2*np.pi #Fuerza de empuje del gas
W_fin=F_a/g #Peso del lecho para caer
#Cambio de masa
#delta_W=W_fin-W
#t_grow=(delta_W)/(r_part*(np.pi*((d_t/2)**2)*L_f))
r_part_reactor=r_part_media*(np.pi*((d_t/2)**2)*L_f) #kg/s de silicio que se depositan en el reactor
#☺print('r_part_reactor (kg/s)', r_part_reactor)
#t_grow=(W_fin-W)/r_part_reactor

#Cantidad del lecho ocupado cada hora 
m_Si=((SiH4_n38/1000)*28.09*X_s*(1-finos_perdidos))
print('m_Si (kg/h)', m_Si)

#Fraccion_lecho=(m_Si*g+W)/W_fin
Fraccion_lecho=(m_Si*g)/(W_fin-W)
print('Fraccion_lecho ', Fraccion_lecho)

m_in=W*Fraccion_lecho
print('m_in (kg/h) ', m_in)





##TIEMPO DE CRECIMIENTO DE UNA PARTICULA
##Masa inicial particulas
#m_ini=(4/3)*np.pi*((d_p/2)**3)*rho_s
##Numero de particulas en el lecho
#n_part=W/(m_ini*5)
##Fuerza de arrastre (depende de Re)
#Re=(d_p*rho_g*u_0)/mi_g
#print ('Re',Re)
##print('Dado que Re>1000 la fuera inercual es la fuera dominante')
##F_a=0.2*rho_g*np.pi*((d_p/2)**2)*(u_0**2) #Componente inercial
##F_a=6*np.pi*mi_g*(d_p/2)*u_0 #Componente de viscosidad
##Empuje
#F_a=delta_P_lecho*(d_t/2)**2*np.pi #Fuerza de empuje del gas
##E=((4/3)*np.pi*((d_p/2)**3))*rho_g*g
##Masa mínima de las partículas para caer
##F_a=F_p=m_fin*g
#m_fin=((F_a)/g)/n_part #Masa de cada particula para caer
#print('m_ini (kg)', m_ini)
#print('m_fin (kg)', m_fin)
##Tiempo de crecimiento
##n_part=(np.pi*((d_t/2)**2)*L_f)/((4/3)*np.pi*((d_p/2)**3))
##t_grow_total=((m_fin-m_ini)*n_part)/(r_part*(np.pi*((d_t/2)**2)*L_f)) #Velocidad de deposicion en todo el reactor
#t_grow=(m_fin-m_ini)/(r_part*(4/3)*np.pi*((d_p/2)**3)) #Velocidad de deposición en una particula
##t_grow=t_grow_total/n_part #Velocidad de deposición en una particula
#print ('Tiempo de crecimiento de las partículas del lecho')
#print ('t_grow (s) = ', t_grow)



##VELOCIDAD DE ENTRADA Y DE SALIDA DE PARTICULAS
##Cada t_grow hay que reemplazar la totalidad de particulas del lecho, por lo que la velocidad de entrada y salida de particulas es la siguiente
#m_in=(W/t_grow)*3600
#m_out=m_in
#print ('Velocidad de entrada y salida de las particulas')
#print ('m_in (kg/h) = ', m_in)
#print ('m_out (kg/h) = ', m_out)
