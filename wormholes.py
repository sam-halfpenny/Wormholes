import math
import numpy
import imageio
import time
from multiprocessing import Pool
import json
from pathlib import Path
FOV=10
#Wormhole centered at origin with throat radius 1
#REMEMBER TO SET LIGHT'S INIT.VEL. TO 1
LD=100000
ds=0.01
univ1=imageio.v2.imread('C:/Users/sam05/OneDrive/Desktop/CODE/github-repositorys/Wormholes/universes/univ1.webp')
univ2=imageio.v2.imread('C:/Users/sam05/OneDrive/Desktop/CODE/github-repositorys/Wormholes/universes/univ2.jpg')
def segScreen(offset):
    workers=offset['workers']
    nim=[]
    for i in range(2*IMSIZE//(workers//2)):
        nim.append([])
        print(str(i))
        for j in range(2*IMSIZE//2):
            raydata=wormholes[0].traceRayFromLookup({'x':-10,'y':0,'z':0},unit_vector({'x':10,'y':-(i/2+offset['x'])*FOV/IMSIZE,'z':(j/2+offset['y'])*FOV/IMSIZE}),1)
            if raydata['univ']==0:
                nim[i].append([0,0,0,255])
            else:
                nim[i].append(findpixel(raydata['vel'],raydata['univ']))
    return nim
def order(arr,workers):
    top=arr[0]
    bottom=arr[1]
    print(len(arr))
    for i in range((workers//2)-1):
        top=numpy.concatenate((top,arr[2*i+2]),axis=0)
        bottom=numpy.concatenate((bottom,arr[2*i+3]),axis=0)
    return numpy.concatenate((bottom,top),axis=1)
def pndiff(p1,p2):
    return {'x':p1['x']-p2['x'],'y':p1['y']-p2['y'],'z':p1['z']-p2['z']}
def unit_vector(a):
    magnitude=math.sqrt(a['x']**2+a['y']**2+a['z']**2)
    if magnitude ==0:
        return {'x':0,'y':0,'z':0}
    a_hat={
        'x':a['x']/magnitude,
        'y':a['y']/magnitude,
        'z':a['z']/magnitude
    }
    return a_hat
def mag(a):
    mag=math.sqrt(a['x']**2+a['y']**2+a['z']**2)
    return mag
def dot_product(v1,v2):
    return v1['x']*v2['x']+v1['y']*v2['y']+v1['z']*v2['z']
def cross_product(a,b):
    final={
        'x':a['y']*b['z']-a['z']*b['y'],
        'y':a['z']*b['x']-a['x']*b['z'],
        'z':a['x']*b['y']-a['y']*b['x']
    }
    return final
def findpixel(vel,univ):
    vel=unit_vector(vel)
    theta=math.atan2(vel['x'],vel['z'])
    phi=math.asin(vel['y'])
    u=(theta/(2*math.pi))+0.5
    v=0.5-(phi/math.pi)
    if univ==1:
        W=len(univ1[0])-1
        H=len(univ1)-1
        pixel=univ1[int(v*H)][int(u*W)].tolist()
        pixel.append(255)
        return pixel
    else:
        W=len(univ2[0])-1
        H=len(univ2)-1
        pixel=univ2[int(v*H)][int(u*W)].tolist()
        pixel.append(255)
        return pixel
class Wormhole:
    def __init__(self,throatRad,pos,traceBoundRad,lookupdensity):
        self.lookupdensity=lookupdensity
        self.throat=throatRad
        self.pos=pos
        self.traceRad=traceBoundRad
        jsonfile='C:/Users/sam05/OneDrive/Desktop/CODE/github-repositorys/Wormholes/lookups/lookup_'+str(self.throat)+'_'+str(self.traceRad)+'_'+str(self.lookupdensity)+'.json'
        my_file = Path(jsonfile)
        if my_file.is_file():
            with open(jsonfile,'r') as jf:
                self.lookup=json.load(jf)
        else:
            self.lookup=self.generateLookup()
    def traceRay(self,pos1,velocity,univ):
        pos1=pndiff(pos1,self.pos)
        vel={'hat':unit_vector(velocity),'mag':mag(velocity),'full':velocity}
        pos={'hat':unit_vector(pos1),'mag':mag(pos1),'full':pos1}
        l_0=math.sqrt(pos['mag']**2-self.throat**2)*univ
        phi_dot=math.sqrt((1-dot_product(vel['hat'],pos['hat'])**2)/((l_0**2+self.throat**2)))
        v_0=numpy.copysign(math.sqrt(1-(l_0**2+self.throat**2)*(phi_dot**2)),dot_product(vel['hat'],pos['hat']))
        b=math.sqrt((l_0**2+self.throat**2)*(1-v_0**2))
        if abs(b-1)<=0.01:
            return {'univ':0}
        phi_0=0
        loop=True
        l=l_0
        v=v_0
        phi=phi_0
        while loop==True:
            v_dot=((b**2)*l)/((l**2+self.throat**2)**2)
            phi_dot=b/(l**2+self.throat**2)
            v=v+ds*v_dot
            l=l+ds*v
            phi=phi+ds*phi_dot
            if math.sqrt(l**2+self.throat**2)>self.traceRad:
                loop=False
        new_mag=math.sqrt(l**2+self.throat**2)
        if dot_product(pos['hat'],vel['hat'])==1:
            new_dir=pos['hat']
            new_vel=new_dir
        else:
            x_dir=pos['hat']
            y_dir=unit_vector({
                'x':vel['hat']['x']-dot_product(vel['hat'],pos['hat'])*pos['hat']['x'],
                'y':vel['hat']['y']-dot_product(vel['hat'],pos['hat'])*pos['hat']['y'],
                'z':vel['hat']['z']-dot_product(vel['hat'],pos['hat'])*pos['hat']['z']
            })
            new_dir={
                'x':math.cos(phi)*x_dir['x']+math.sin(phi)*y_dir['x'],
                'y':math.cos(phi)*x_dir['y']+math.sin(phi)*y_dir['y'],
                'z':math.cos(phi)*x_dir['z']+math.sin(phi)*y_dir['z']
            }
            new_tan={
                'x':math.cos(phi)*y_dir['x']-math.sin(phi)*x_dir['x'],
                'y':math.cos(phi)*y_dir['y']-math.sin(phi)*x_dir['y'],
                'z':math.cos(phi)*y_dir['z']-math.sin(phi)*x_dir['z']
            }
            r_dot=(v*l)/math.sqrt(l**2+self.throat**2)
            new_vel={
                'x':r_dot*new_dir['x']+math.sqrt(1-r_dot**2)*new_tan['x'],
                'y':r_dot*new_dir['y']+math.sqrt(1-r_dot**2)*new_tan['y'],
                'z':r_dot*new_dir['z']+math.sqrt(1-r_dot**2)*new_tan['z'],
            }
        univ=l/abs(l)
        new_pos={
            'x':self.pos['x']+new_dir['x']*new_mag,
            'y':self.pos['y']+new_dir['y']*new_mag,
            'z':self.pos['z']+new_dir['z']*new_mag
        }
        return {'pos':new_pos,'vel':new_vel,'univ':univ}
    def traceRayFromLookup(self,pos1,velocity,univ):
        pos1=pndiff(pos1,self.pos)
        vel={'hat':unit_vector(velocity),'mag':mag(velocity),'full':velocity}
        pos={'hat':unit_vector(pos1),'mag':mag(pos1),'full':pos1}
        if abs(pos['mag']-self.traceRad)>0.05:
            print('E: cannot lookup non-surface paths')
            return 'ERROR'
        theta=math.acos(-dot_product(vel['hat'],pos['hat']))
        i=round((2*theta*self.lookupdensity)/math.pi)
        data=self.lookup[i][1]
        if data['univ']==0:
            return {'univ':0}
        phi=data['phi']
        l=data['l']
        v=data['v']

        new_mag=math.sqrt(l**2+self.throat**2)
        if dot_product(pos['hat'],vel['hat'])==1:
            new_dir=pos['hat']
            new_vel=new_dir
        else:
            x_dir=pos['hat']
            y_dir=unit_vector({
                'x':vel['hat']['x']-dot_product(vel['hat'],pos['hat'])*pos['hat']['x'],
                'y':vel['hat']['y']-dot_product(vel['hat'],pos['hat'])*pos['hat']['y'],
                'z':vel['hat']['z']-dot_product(vel['hat'],pos['hat'])*pos['hat']['z']
            })
            new_dir={
                'x':math.cos(phi)*x_dir['x']+math.sin(phi)*y_dir['x'],
                'y':math.cos(phi)*x_dir['y']+math.sin(phi)*y_dir['y'],
                'z':math.cos(phi)*x_dir['z']+math.sin(phi)*y_dir['z']
            }
            new_tan={
                'x':math.cos(phi)*y_dir['x']-math.sin(phi)*x_dir['x'],
                'y':math.cos(phi)*y_dir['y']-math.sin(phi)*x_dir['y'],
                'z':math.cos(phi)*y_dir['z']-math.sin(phi)*x_dir['z']
            }
            r_dot=(v*l)/math.sqrt(l**2+self.throat**2)
            new_vel={
                'x':r_dot*new_dir['x']+math.sqrt(1-r_dot**2)*new_tan['x'],
                'y':r_dot*new_dir['y']+math.sqrt(1-r_dot**2)*new_tan['y'],
                'z':r_dot*new_dir['z']+math.sqrt(1-r_dot**2)*new_tan['z'],
            }
        new_pos={
            'x':self.pos['x']+new_dir['x']*new_mag,
            'y':self.pos['y']+new_dir['y']*new_mag,
            'z':self.pos['z']+new_dir['z']*new_mag
        }
        return {'pos':new_pos,'vel':new_vel,'univ':univ*data['univ']}
    def generateLookup(self):
        jsonoutfile='C:/Users/sam05/OneDrive/Desktop/CODE/github-repositorys/Wormholes/lookups/lookup_'+str(self.throat)+'_'+str(self.traceRad)+'_'+str(self.lookupdensity)+'.json'
        jof=open(jsonoutfile,'w')
        data=[]
        for t in range(self.lookupdensity):
            print(t)
            theta=(t*math.pi)/(self.lookupdensity*2)
            l_0=math.sqrt((self.traceRad-0.001)**2-self.throat**2)
            phi_dot=math.sqrt((1-math.cos(theta)**2)/((l_0**2+self.throat**2)))
            v_0=-math.sqrt(1-(l_0**2+self.throat**2)*(phi_dot**2))
            b=math.sqrt((l_0**2+self.throat**2)*(1-v_0**2))
            if abs(b-1)<=0.0001:
                data.append([t,{'univ':0}])
            else:
                phi_0=0
                loop=True
                l=l_0
                v=v_0
                phi=phi_0
                count=0
                while loop==True:
                    count+=1
                    v_dot=((b**2)*l)/((l**2+self.throat**2)**2)
                    phi_dot=b/(l**2+self.throat**2)
                    v=v+ds*v_dot
                    l=l+ds*v
                    phi=phi+ds*phi_dot
                    if math.sqrt(l**2+self.throat**2)>self.traceRad:
                        loop=False
                data.append([t,{'phi':phi,'univ':abs(l)/l,'v':v,'l':l}])
        jof.write(json.dumps(data))
        return data
    


nimpath='C:/Users/sam05/OneDrive/Desktop/CODE/github-repositorys/Wormholes/render/Wormhole'+str(round(time.time()))+'.png'
wormholes=[Wormhole(1,{'x':0,'y':0,'z':0},10,LD)]
pos={'x':5,'y':1,'z':0}
vel=unit_vector({'x':-7,'y':1,'z':0})
IMSIZE=2000
nim=[]
worker_count=16
args=[]
if __name__ == '__main__':
    # arr=multiprocessing.Array('i',range(workerpool))
    for i in range(worker_count):
        args.append({'x':(i//2)*(IMSIZE/(worker_count/2))-IMSIZE/2,'y':(i%2)*(-IMSIZE/2),'workers':worker_count})
    with Pool(worker_count) as p:
        print(worker_count,args)
        nim=order(p.map(segScreen, args),worker_count)
    print(nim)
    nim = numpy.clip(nim, 0, 255).astype(numpy.uint8)
    imageio.imwrite(nimpath,nim)
# for i in range(IMSIZE):
#     print(i)
#     nim.append([])
#     for j in range(IMSIZE):
#         raydata=wormholes[0].traceRay({'x':0,'y':0,'z':-8},unit_vector({'x':i*10/IMSIZE,'y':j*10/IMSIZE,'z':10}),1)
#         if raydata['univ']==0:
#             nim[i].append([255,0,0,255])
#         else:
#             nim[i].append(findpixel(raydata['vel'],raydata['univ']))





