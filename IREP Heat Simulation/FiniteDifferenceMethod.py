import numpy as np
import pygame
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import asyncio
import matplotlib.backends.backend_agg as agg
import pylab
from pygame import * 
from random import randint
from math import lcm
from matplotlib.animation import FuncAnimation 
from matplotlib.pyplot import imshow
from matplotlib import cm
from decimal import *
from pathlib import Path

#Initiation Procedures
pygame.init()
clock = pygame.time.Clock()
pygame.font.init() 

#Generation File 
GenerationFile = [
    {
        'Session Name': 'CopperTest',
        'height': 42,
        'unitsize': 6,
        #Boudary Conditions
        'top': 200, 
        'left': 25, 
        'right': 25, 
        'down': 25, 
        #Iteration Time 
        'iterationtime':0.9,    
        #Snpashot:
        'frequency': 0.1, 
        'snapshotpos': 10, 
        'maxgraphtemp': 220, 
        'mingraphtemp': 20
    },
    
    {
        'name': 'Copper', 
        'property': 111, 
        'length': 32
    },



]



class Simulation(): 
    def __init__(self, genfile):
    #Loading Frames
        self.frames = []
        self.STEP = 1 



    ## General Settings
        self.genfile = genfile
        generalsettings = genfile[0]
        self.ems = [genfile[x]['property'] for x in range(1, len(genfile))]
    
        #Size of the simulation
        self.unitsize = generalsettings['unitsize']
         
        #Boundary Conditions
        self.topbound = generalsettings['top']
        self.bottombound = generalsettings['down']
        self.leftbound = generalsettings['left']
        self.rightbound = generalsettings['right']

        self.bounds = max([self.topbound, self.bottombound, self.leftbound, self.rightbound])

        #Time/Space Conditions

        self.deltax = 1
        self.elapsed = 0
        self.iterationtime = generalsettings['iterationtime']

        #Snapshot 

        self.snapshotpos = generalsettings['snapshotpos'] 
        self.maxgraphtemp = generalsettings['maxgraphtemp']
        self.mingraphtemp = generalsettings['mingraphtemp']
        self.frequency = generalsettings['frequency']

    ## Layers and Dimensions
        self.PADDING = 50 

        self.height = generalsettings['height']
        self.length = sum([genfile[x]['length'] for x in range(1, len(genfile))])



    ##Setting up the array 

        #Surface Array
        self.Surface = Surface((self.length * self.unitsize + self.PADDING * 2, self.height * self.unitsize + self.PADDING))
        self.simSurface = Surface((self.length * self.unitsize, self.height * self.unitsize))
        self.Surface.fill(WHITE)
        self.simSurface.fill(WHITE)
        self.surface = np.empty((self.height, self.length))
        self.surface.fill(25)
        self.surface[:1, :,] = self.topbound
        self.surface[self.height - 1:, :,] = self.bottombound
        self.surface[:, :1,] = self.leftbound
        self.surface[:, self.length - 1:,] = self.rightbound

        #Alpha Array 
        alphas = [] 
        positionlist = [0]

    
        for i in range(1, len(genfile)): 
            alpha = np.empty((self.height, int(genfile[i]['length'])))
            alpha.fill(genfile[i]['property'])
            alphas.append(alpha)
        
        self.alpha = np.concatenate(tuple(alphas), axis = 1)

        #Delta T Array 

        self.deltat = (self.deltax ** 2)/ (4 * self.alpha)
        rawframecycle = list(dict.fromkeys(list(self.deltat[0])))
        self.rawframecycle = [float(rawframecycle[x]) for x in range(len(rawframecycle))]

        self.rawframecycle.sort()
    
        minvalue = min(rawframecycle)
        i = 0 
        while minvalue < 1: 
            i += 1
            minvalue = minvalue * 10 
        
        self.roundedt = np.around(self.deltat * (10 ** i))
        
        #This part assumes that we are only concerned with simulations that are layered horizontally 

        roundframecycle = list(dict.fromkeys(list(self.roundedt[0])))
        self.framecycle = [int(roundframecycle[x]) for x in range(len(roundframecycle))]
        self.framecycle.sort()
        self.minStep = min(self.framecycle)
        self.minDeltaT = min(self.rawframecycle)
        
        #framerate = np.lcm.reduce(framecycle)
        self.linePos = [0]
        for i in range(1, len(GenerationFile)): 
            newPos = GenerationFile[i]['length'] + self.linePos[i-1]
            self.linePos.append(newPos)


    def step(self): 
        for i in range(1, self.height -1, self.deltax): 
            for j in range(1, self.length -1, self.deltax):
                if self.STEP % self.roundedt[i][j] == 0: 
                    alpha = self.alpha[i][j]
                    deltat = self.deltat[i][j]

                    condition = (alpha * deltat) / (self.deltax ** 2)
                    self.surface[i, j] = condition * (self.surface[i+1][j] + self.surface[i-1][j] + self.surface[i][j+1] + self.surface[i][j-1] - 4*self.surface[i][j]) + self.surface[i][j]
        self.STEP += 1
         
        if self.STEP % self.minStep == 0: 
            self.elapsed += self.minDeltaT

        

    def draw(self): 
        self.units = []
        self.lines = []
        self.markers = []
        self.textmarkers = []
        
        font = pygame.font.Font('cour.ttf', 20)

        for i in range(len(self.surface)): 
            for j in range(len(self.surface[i])): 
                self.units.append(
                    Unit(
                        state = self.surface[i][j], 
                        bound = 200, 
                        coord = (j * self.unitsize, i * self.unitsize), 
                        unitsize = self.unitsize 

                    )
                )
    
        for unit in self.units: 
            unit.draw(self.simSurface)
    
        for ems in self.ems: 
            layerData = Data(
                fontsize = 12,
                data = [
                    f'ε = {ems}'
                ]
            ) 
            textmarker = layerData.draw()
            self.textmarkers.append(textmarker)



        for i in range(len(self.linePos)): 
            overlay = Overlay(
                start = (self.linePos[i] * self.unitsize, 0),
                end = (self.linePos[i] * self.unitsize, self.height * self.unitsize)
            )
        
                
            self.lines.append(overlay)
            self.markers.append(font.render(str(self.linePos[i]), True, BLACK))        

        for i in range(len(self.lines)): 
            self.lines[i].draw(self.simSurface)

        heightmarker = font.render(str(self.height), True, BLACK)
        snapshotmarker = font.render('<-', True, BLACK)

        for i in range(len(self.markers)): 
            markerposition = self.linePos[i] * self.unitsize + self.PADDING / 2 
            self.Surface.blit(self.markers[i], (markerposition, self.height * self.unitsize + self.PADDING/2))
        for i in range(len(self.markers) - 1):
            markerposition = self.linePos[i] * self.unitsize + self.PADDING / 2 
            self.simSurface.blit(self.textmarkers[i], (markerposition - (4 * self.unitsize), (self.height - 3) * self.unitsize))

        self.Surface.blit(self.simSurface, (self.PADDING / 2, self.PADDING / 2))
        self.Surface.blit(heightmarker, (0, self.PADDING/2))
        self.Surface.blit(snapshotmarker, (self.simSurface.get_width() + self.PADDING/2, self.PADDING/2 + self.unitsize * (self.snapshotpos - 2)))

        

        return self.Surface 


    def snapshot(self): 
        f= open('snapshotdata.txt', 'a+')

        self.position = self.snapshotpos # At which point is the snapshot taken 
        snapshot = list(self.surface[self.position])
        self.x = [x for x in range(len(snapshot))]

        for i in range(len(self.linePos) - 1): 
            yaxis = 140
            arrowstart = self.linePos[i]
            arrowend = self.linePos[i + 1]
            plt.arrow(arrowstart, yaxis, arrowend - arrowstart, 0, head_width=2, head_length=0.5, linewidth=0.75, color='black', length_includes_head=True)
            plt.arrow(arrowend, yaxis, arrowstart - arrowend, 0, head_width=2, head_length=0.5, linewidth=0.75, color='black', length_includes_head=True)


        for line in self.linePos: 
            plt.axvline(x = line, color = 'black', label = 'axvline - full height', linestyle = '--', ymax = 1, linewidth = 0.5)
        
        plt.plot(self.x, snapshot, color = 'blue', linewidth = 0.5)
        plt.xlim([0, self.length])
        plt.ylim ([self.mingraphtemp, self.maxgraphtemp])
        plt.xlabel("Material (mm)", family='Courier New', size=14)
        plt.ylabel("Temperature (C)", family='Courier New', size=14)
        plt.savefig('snapshot.png', pad_inches = 0, bbox_inches = 'tight')
        for i in range(len(self.x)): 
            f.write(f'{self.x[i]} {snapshot[i]} \n')

        f.write('\n')
        f.close()

class Unit(): 
    def __init__(self, coord, unitsize, state, bound):
        self.bound = bound
        self.state = state
        self.coord = coord  
        self.x = self.coord[0]
        self.y = self.coord[1]
        self.unitsize = unitsize
        

    def draw(self, win): 
        colortuple = cm.jet(self.state / self.bound)
        self.color = tuple([colortuple[x] * 255 for x in range(len(colortuple))])
        
        # self.color = colorspectrum(self.state, self.bound)
        pygame.draw.rect(
            surface = win, 
            color = self.color, 
            rect = (self.x, self.y, self.unitsize, self.unitsize),
        )

class Overlay():
    def __init__(self, start, end): 
        self.start = start, 
        self.end = end
        
    def draw(self, surface): 
        self.surface = surface
        pygame.draw.line(
            surface = self.surface, 
            start_pos = self.start, 
            end_pos = self.end, 
            color = BLACK
        )
        
class Data(): 
    def __init__(self, data, fontsize):
        self.fontsize = fontsize  , 
        self.data = data
        font = pygame.font.Font('cour.ttf', fontsize)
        self.text = [font.render(self.data[x], True, BLACK) for x in range(len(self.data))]
        lenchar = [len(self.data[x]) for x in range(len(self.data))]
        maxchar = lenchar.index(max(lenchar))

        self.charwidth = self.text[maxchar].get_width() / max(lenchar)
        self.charheight = self.text[0].get_height()


        self.length = self.text[maxchar].get_width() + 2*self.charwidth
        self.width = self.charheight * len(self.data)
        self.Surface = Surface((self.length, self.width))
        self.Surface.fill(WHITE)
        pygame.draw.rect(
            surface = self.Surface,
            color = BLACK, 
            rect = (0, 0, self.length, self.width),
            width = 1,
            
        )

        
    
    def draw(self):
        for i in range(len(self.text)):
            self.Surface.blit(self.text[i], (self.charwidth, i * self.charheight))

        return self.Surface 
#Color Constants

GREY = (100, 100, 100)
BLACK = (0, 0, 0)
WHITE = (255, 255, 255)


def main(): 
    #Simulation initialized
    sim = Simulation(
            genfile = GenerationFile, 

        )

    sessionname = sim.genfile[0]['Session Name']
    simheight = sim.draw().get_height()
    simwidth = sim.draw().get_width()
    iterationtime = sim.iterationtime

    #Initial Figure Setup 
    f= open('snapshotdata.txt', 'w+')
    f.write(' ')
    f.close()
    sim.snapshot()
    image = pygame.image.load('snapshot.png')
    

    plt.figure()
    a = np.array([[0, sim.bounds]])
    pylab.figure(figsize=(simwidth / 98, 0.5))
    img = pylab.imshow(a, cmap="jet")
    pylab.gca().set_visible(False)
    cax = pylab.axes([0.1, 0.2, 0.8, 0.6])
    pylab.colorbar(orientation="horizontal", cax=cax)
    pylab.savefig("colorbar.png", bbox_inches='tight', pad_inches = 0)
    plt.figure()

    colorbar = pygame.image.load('colorbar.png')

    #Initial Text display Setup

    layers = []
    for i in range(1, len(GenerationFile)): 
        name = GenerationFile[i]['name']
        property = GenerationFile[i]['property']
        length = GenerationFile[i]['length']
        layer = [
            f'Layer: {name}', 
            f'Diffusivity: {property}', 
            f'Length: {length}', 
            f'                      '
            
            ]
        for ele in layer: 
            layers.append(ele)
        

    layersdata = Data(
        fontsize = 18, 
        data = [
            f'Layers: {len(GenerationFile) - 1}', 
            f'---------------------------------'                               
        ] + layers
    )

    snapfreq = GenerationFile[0]['frequency']
    graphdata = Data(
            fontsize = 18, 
            data = [
                'Graph Settings',
                f'-----------------------------------', 
                f'Snapshot Position (horizontal): {sim.snapshotpos}',
                f'Maximum Graph Temperature: {sim.maxgraphtemp}',
                f'Minimum Graph Temperature: {sim.mingraphtemp}',
                f'Snapshot Frequency: {snapfreq} s',
            ]
        )

    image = pygame.image.load('snapshot.png')
    imagesurface = Surface((image.get_width(), image.get_height()))
    imagesurface.blit(image, (0, 0))
    figurenote = Data(
        fontsize = 12, 
        data = [
            f'Duration = {sim.iterationtime}s',
            f'Freq = {sim.frequency}s', 
            f'Position = {sim.snapshotpos} mm'
        ]
    )

    timeandBoundary = Data(
            fontsize = 18,
            data = [ 
                f'Simulation Duration: {iterationtime} s', 
                'Boundary Conditions: ', 
                f'   Top: {sim.topbound}', 
                f'   Bottom: {sim.bottombound}', 
                f'   Left: {sim.leftbound}', 
                f'   Right: {sim.rightbound}', 
                f'                                       ',
                f'---------------------------------------'
            ]
        )
    imagesurface.blit(figurenote.draw(), (imagesurface.get_width() * 0.15, 5))

    #Snapshot Initialization
    snapshots = []
    times = []
    decfreq = Decimal(str(sim.frequency))
    decitertime = Decimal(str(sim.iterationtime))

    for i in range(int(decitertime / decfreq) + 1): 
        snapshots.append(float(str(i * decfreq)))
    # Window Settings

    defaultlength = 1300
    LENSCALE = (simwidth + image.get_width())/1600 
    HEIGHTSCALE = (image.get_height() + layersdata.width + 20)/900

    if 1600 * LENSCALE + 100 < defaultlength:
        WINDIMENSION = (defaultlength, 900 * HEIGHTSCALE + 50)
    else: 
        WINDIMENSION = (1600 * LENSCALE + 100, 900 * HEIGHTSCALE + 50)

    WINDOW = pygame.display.set_mode(WINDIMENSION,HWSURFACE|DOUBLEBUF|RESIZABLE)
    
    DISPLAY = WINDOW.copy()
    pygame.display.set_caption('Finite Difference Method Simulator')
    DISPLAY.fill(WHITE)

    time = 0
    SIMSURFACE = Surface((simwidth, simheight + colorbar.get_height() + 40))  
    SIMSURFACE.fill(WHITE) 


    #Simulation Loop
    while time < iterationtime: 
        
        for event in pygame.event.get():
            if event.type == pygame.QUIT:
                running = False
                pygame.QUIT()

        #Timebox generation
        stringtime = '{:.3f}'.format(round(sim.elapsed, 5))
        time = float(stringtime)
        stringtime = f'{stringtime} s' 

        timebox = Data(
            fontsize = 18,
            data = [
                f'Elapsed Time: {stringtime}'
            ]
        )

        #UI Blitting

        SIMSURFACE.blit(sim.draw(), (0, 0))
        SIMSURFACE.blit(timebox.draw(), (sim.PADDING / 2, simheight + colorbar.get_height() + 10))
        SIMSURFACE.blit(colorbar, (sim.PADDING / 2, simheight ))

        DISPLAY.blit(SIMSURFACE, (0, 0))
        DISPLAY.blit(timeandBoundary.draw(), (sim.PADDING / 2, imagesurface.get_height() + 10))
        DISPLAY.blit(layersdata.draw(), (sim.PADDING / 2 + timeandBoundary.length + 10, imagesurface.get_height() + 10))
        DISPLAY.blit(graphdata.draw(), (sim.PADDING / 2 + timeandBoundary.length + layersdata.length + 20, imagesurface.get_height() + 10))
        DISPLAY.blit(imagesurface, (simwidth + 20, 0))

        #Snapshot Conditions & Simulation progression

        if time <= iterationtime:
            if len(snapshots) > 0:
                nextshot = snapshots[0]
                if time == nextshot and time not in times:
                    pygame.image.save(SIMSURFACE, f'C:/Users/rythc/VSCode Python/IREP Heat Simulation/ScreenShots/{sessionname}@{time}s.png')
                    sim.snapshot()
                    snapshots.pop(0)
                    image = pygame.image.load('snapshot.png')
                    imagesurface = Surface((image.get_width(), image.get_height()))
                    imagesurface.blit(image, (0, 0))
                    DISPLAY.blit(imagesurface, (simwidth + 20, 0))
                    sim.step()
                elif time + sim.minDeltaT > nextshot and time not in times:
                    pygame.image.save(SIMSURFACE, f'C:/Users/rythc/VSCode Python/IREP Heat Simulation/ScreenShots/{sessionname}@{time}s.png')
                    sim.snapshot()
                    snapshots.pop(0)
                    image = pygame.image.load('snapshot.png')
                    imagesurface = Surface((image.get_width(), image.get_height()))
                    imagesurface.blit(image, (0, 0))
                    DISPLAY.blit(imagesurface, (simwidth + 20, 0))
                    sim.step()
                else: 
                    sim.step()
            else: 
                sim.step()

        #Save Conditions

        
        

        if time == iterationtime: 
            pygame.image.save(imagesurface, f'C:/Users/rythc/VSCode Python/IREP Heat Simulation/ScreenShots/Figure.jpg')
            converged = (sim.surface[sim.position])
            f= open('converged.txt', 'w+')
            for i in range(len(converged)): 
                f.write(f'{i} {converged[i]}\n')
                
            f.close()

            plt.figure()
            plt.plot([x for x in range(len(converged))], converged)
            plt.savefig('converged.png')
            
            

            break
       
        #Frame Elapsing

        times.append(time)
        times = list(dict.fromkeys(times))
        WINDOW.blit(pygame.transform.scale(DISPLAY, WINDOW.get_rect().size), (0, 0))
        clock.tick(60)
        pygame.display.flip()

    
    #Run-off code

    while True: 
        for event in pygame.event.get():
            if event.type == pygame.QUIT:
                running = False
                pygame.QUIT()

        WINDOW.blit(pygame.transform.scale(DISPLAY, WINDOW.get_rect().size), (0, 0))
        clock.tick(60)
        pygame.display.flip()

    
main()



