import sys
import h5py
#import os
import numpy as np
from random import randint

from bokeh.models.widgets import CheckboxButtonGroup, TextInput, Select, RadioButtonGroup
from bokeh.io import show, curdoc
from bokeh.plotting import figure
from bokeh.layouts import layout, widgetbox, row, gridplot, column
from bokeh.models import LinearAxis, Range1d, HoverTool, ColumnDataSource

INTERVAL = 1000
SPLIT = " "
filesList = []

#Einlesen der ConfigDatei und Abspeichern der einzelnen Dateien in filesList
if (len(sys.argv) == 2):
    try:
        #TODO: Prüfen ob PFad existiert? pathExists soll es als Methode gegeben
        configFile = open(sys.argv[1], 'r')
        configFileContent = configFile.readlines()
        #lesen der Zeilen; startet bei der ersten Zeile, den Überschriften
        for x in range (len(configFileContent)):
            #Sample, Class, Raw, Derived: 4 Spalten
            #Spalten bei dem gegebenen SPLIT Punkt
            line = configFileContent[x].split(SPLIT)
            #in FilesList eintragen
            filesList.append(line)
    except:
        print("Etwas ist schief gelaufen, wahrscheinlich hat die config-Datei nicht das richtige Format")
else:
    print("Bitte die korrekte Anzahl an Argumenten übergeben: bokeh serve --show {Dateinamen} --args {Pfad zur Config Datei}")

#Definieren der indexFile
indexFile = h5py.File(filesList[1][2], 'r')

#Berechnen der Level pro Sample, wenn Index und Daten gegeben sind
def calculateMethylationLevel(index, data, sampleNameList):
    #Berechnung der einzelnen Rohdaten
    dataList = []
    dataValues={'sampleName': [], 'index': [], 'methylatedLevel': [], 'methylatedAbsolut' : [], 'errorValue' : [], 'methylationErrorStart' : [], 'methylationErrorEnd' : []}
    #Fasse alle Daten anhand der gewuenschten Kategorie (Sample, Klasse) zusammen: Bei gleichem Index müssen methyliertes bzw. unmethyliertes zusammengerechnet werden
    dataValueSumed = {'sampleName': [], 'index':[], 'methylated': [], 'unmethylated': []}
    dataSumed = []
    for x in range(0, len(data)):
        sample = sampleNameList[x]
        found = False
        foundIndex = -1
        methylationMethylatedSum = []
        methylationUnmethylatedSum = []
        sampleName = []
        #Sucht, ob entsprechendes Sample/Klasse schon vorhanden ist
        for m in range(0, len(dataSumed)):
            if sample in dataSumed[m]['sampleName']:
                found = True
                foundIndex = m
        for y in range(0, len(data[x])):
            # finden der gleichen Indices und addieren der Werte
            if(found):
                previousValue = dataSumed[foundIndex]['methylated'][y]
                nextValue = data[x][y][0]+data[x][y][2]
                newValue = previousValue+nextValue
                dataSumed[foundIndex]['methylated'][y] = newValue
                previousValue = dataSumed[foundIndex]['unmethylated'][y]
                nextValue = data[x][y][1]+data[x][y][3]
                newValue = previousValue+nextValue
                dataSumed[foundIndex]['unmethylated'][y] = newValue
            # sonst wurde ein neues Sample/Klasse gefunden und muss addiert werden
            else:
                sampleName.append(sample)
                methylated = data[x][y][0]+data[x][y][2]
                unmethylated = data[x][y][1]+data[x][y][3]
                methylationMethylatedSum.append(methylated)
                methylationUnmethylatedSum.append(unmethylated)
        #Datensatz muss noch eingefuegt werden, wenn nicht gefunden
        if(not found):
            dataValues={'sampleName': sampleName, 'index': index, 'methylated':methylationMethylatedSum, 'unmethylated':methylationUnmethylatedSum}
            dataSumed.append(dataValues)
    #Neuberechnung der Methylierungslevel
    # Gehe alle Eintraege der zusammengefassen Daten durch
    for x in range(0, len(dataSumed)):
        methylationSum = []
        error = []
        startError = []
        endError = []
        methylationLevel = []
        methylationErrorStart = []
        methylationErrorEnd = []
        # Berechnung fuer die einzelnen Datenpunkte
        for y in range(0, len(dataSumed[x]['index'])):
            methylationSum = dataSumed[x]['methylated'][y]+dataSumed[x]['unmethylated'][y]
            if (methylationSum == 0):
                methylationLevel.append(0)
                startError.append(0)
                endError.append(0)
                error.append(0)
            else:
                level = dataSumed[x]['methylated'][y]/methylationSum
                sigma = ((level)*(1-level))/np.sqrt(methylationSum)
                error.append(sigma)
                startError.append((level-sigma))
                endError.append((level+sigma))
                methylationLevel.append(level)
        dataValues = {'sampleName': dataSumed[x]['sampleName'], 'index': index, 'methylationLevel': methylationLevel, 'methylatedAbsolut' : dataSumed[x]['methylated'], 'errorValue' : error, 'methylationErrorStart' : startError, 'methylationErrorEnd' : endError}
        dataList.append(dataValues)
    return dataList

#Pre: chromosom ist im Stringformat, start und end sind im Integerformat, sample ist ein Boolean und gibt an, ob pro Sample (true)
# oder pro Klasse (false) angezeigt werden soll, difference ist ein Boolean und gibt an, ob der Unterschied zwischen den Klassen hervorgehoben werden soll
def displayDataRaw(chromosom, start, end, sample, difference):
    global filesList, finalGridplot, p
    #Ueberpruefe, ob das Intervall in den gegeben Grenzen liegt
    namesSampleClass = []
    if (start >= 0 and end > 0 and start <=end and end-start<= INTERVAL):
        #Hole Daten von allen Samples, ersten beiden Eintraege muessen ignoriert werden
        data = []
        for i in filesList[2:]:
            #Einlesen und Anhaengen der Methylierungen
            data.append(h5py.File(i[2], 'r').get(chromosom).get('methylation')[start:end+1])
            # wenn Sample aktiv ist, einfach nur pro Sample einlesen
            if(sample):
                namesSampleClass.append(i[0])
            else:
                namesSampleClass.append(i[1])
        index = list(indexFile.get(chromosom).get('cpg_positions')[start:end+1])
    else:
        print("Something went wrong")
    # Data hat das folgende Format: pro Sample/Klasse{'sampleName', 'index', 'methylationLevel', 'methylatedAbsolut', 'errorValue' , 'methylationErrorStart', 'methylationErrorEnd' }
    dataCalculated = calculateMethylationLevel(index, data, namesSampleClass)
    #Generate Colors
    # TODO diese statisch machen, sodass es nicht immer wechselt
    color = []
    for i in range(len(dataCalculated)):
        colorRandom = (randint(0,255), randint(0,255), randint(0,255))
        color.append(colorRandom)
    # Gehe Samples durch
    methylationLevelAll = []
    # Alle Sample oder Klassen durchgehen und Punkte zeichnen
    for i in range (len(dataCalculated)):
        source = ColumnDataSource(data=dataCalculated[i])
        colorCircle = (randint(0,255), randint(0,255), randint(0,255))
        colorBar = (randint(0,255), randint(0,255), randint(0,255))
        points = p.circle('index', 'methylationLevel', source=source, size=10, color=color[i], legend=dataCalculated[i]['sampleName'][0])
        error = p.vbar(x='index', width=0.5, bottom='methylationErrorStart', top='methylationErrorEnd', color="black", source=source)
        bars = p.vbar(x='index', width=5, bottom=0, top='methylatedAbsolut', color=color[i], y_range_name="absoluteData", fill_alpha=0.10, source=source)
        p.add_tools(HoverTool(renderers=[points], tooltips=[("index","@index"), ("methylationLevel","@methylationLevel"), ("sample", "@sampleName")]))
        p.add_tools(HoverTool(renderers=[bars], tooltips=[("index", "@index"), ("methylationAbsolut", "@methylatedAbsolut"), ("sample", "@sampleName")]))
        methylationLevelAll.append(dataCalculated[i]['methylationLevel'])
    #Wenn Difference gefordert (Dazu sample auf false, damit nur KLassenunterschiede angezeigt werden), muss diese noch ermittelt werden
    if(difference):
        methylationMin = []
        methylationMax = []
        indexReverse = []
        minReverse = []
        # alle Datenpunkte herausfinden
        for x in range(len(methylationLevelAll[0])):
            # pro Datensatz an der Stelle x herausfinden, was das Minimum und was das Maximum ist
            min = 2
            max = -1
            # alle Datensätze an der Stelle durchgehen
            for y in range (len(methylationLevelAll)):
                if methylationLevelAll[y][x] < min:
                    min = methylationLevelAll[y][x]
                if methylationLevelAll[y][x] > max:
                    max = methylationLevelAll[y][x]
            methylationMin.append(min)
            methylationMax.append(max)
        indexReverse = dataCalculated[0]['index'][:]
        indexReverse.reverse()
        indexOrdered = dataCalculated[0]['index']+indexReverse
        minReverse = methylationMin
        minReverse.reverse()
        methylationMix = methylationMax+ minReverse
        p.patch(indexOrdered, methylationMix, alpha=0.5, line_width=2)

    #Berechnung fuer das Maximum der absoluten Werte fuer eine bessere y-Achsen Darstellung
    maxAbsolut=0
    for x in range (0, len(dataCalculated)):
        maximumYRange = list(dataCalculated[x]['methylatedAbsolut'])
        if (np.max(maximumYRange) > maxAbsolut):
            maxAbsolut = np.max(maximumYRange)
    p.extra_y_ranges ={"absoluteData":Range1d(start=0, end=2*maxAbsolut)} #max(data[x].methylatedAbsolut)+1
    p.y_range=Range1d(0,1.1)
    p.add_layout(LinearAxis(y_range_name="absoluteData"), 'right')

#Pre: chromosom ist im Stringformat, start und end sind im Integerformat, mode hat entweder den Wert 2 für Rohdaten oder 3 für Klassen
def displayDataAlgorithm(chromosom, start, end, sample, difference):
    global filesList, finalGridplot, p
    #Ueberpruefe, ob das Intervall in den gegeben Grenzen liegt
    sampleNames = []
    if (start >= 0 and end > 0 and start <=end and end-start<= INTERVAL):
        #Hole Daten von allen Samples, ersten beiden Eintraege muessen ignoriert werden
        data = []
        for i in filesList[2:]:
            #Einlesen und Anhaengen der Methylierungen
            data.append(h5py.File(i[3], 'r').get(chromosom).get('methylation')[start:end+1])
            # wenn Sample aktiv ist, einfach nur pro Sample einlesen
            if(sample):
                sampleNames.append(i[0])
            else:
                sampleNames.append(i[1])
        index = list(indexFile.get(chromosom).get('cpg_positions')[start:end+1])
    else:
        print("Something went wrong")
    # Data hat das folgende Format: pro Sample{'sampleName', 'index', 'methylationLevel', 'methylatedAbsolut', 'errorValue' , 'methylationErrorStart', 'methylationErrorEnd' }
    dataCalculated = calculateMethylationLevel(index, data, sampleNames)
    #Generate Colors
    # TODO diese statisch machen, sodass es nicht immer wechselt
    color = []
    for i in range(len(dataCalculated)):
        colorRandom = (randint(0,255), randint(0,255), randint(0,255))
        color.append(colorRandom)
    # Gehe Samples durch
    methylationLevelAll = []
    # Alle Sample oder Klassen durchgehen und Punkte zeichnen
    for i in range (len(dataCalculated)):
        source = ColumnDataSource(data=dataCalculated[i])
        colorCircle = (randint(0,255), randint(0,255), randint(0,255))
        colorBar = (randint(0,255), randint(0,255), randint(0,255))
        points = p.circle('index', 'methylationLevel', source=source, size=10, color=color[i], legend="Algorithmus " + dataCalculated[i]['sampleName'][0])
        p.add_tools(HoverTool(renderers=[points], tooltips=[("index","@index"), ("methylationLevel","@methylationLevel"), ("sample", "Alg "+"@sampleName")]))
        methylationLevelAll.append(dataCalculated[i]['methylationLevel'])
    #Wenn Difference gefordert (Dazu sample auf false, damit nur KLassenunterschiede angezeigt werden), muss diese noch ermittelt werden
    if(difference):
        methylationMin = []
        methylationMax = []
        indexReverse = []
        minReverse = []
        # alle Datenpunkte herausfinden
        for x in range(len(methylationLevelAll[0])):
            # pro Datensatz an der Stelle x herausfinden, was das Minimum und was das Maximum ist
            min = 2
            max = -1
            # alle Datensätze an der Stelle durchgehen
            for y in range (len(methylationLevelAll)):
                if methylationLevelAll[y][x] < min:
                    min = methylationLevelAll[y][x]
                if methylationLevelAll[y][x] > max:
                    max = methylationLevelAll[y][x]
            methylationMin.append(min)
            methylationMax.append(max)
        indexReverse = dataCalculated[0]['index'][:]
        indexReverse.reverse()
        indexOrdered = dataCalculated[0]['index']+indexReverse
        minReverse = methylationMin
        minReverse.reverse()
        methylationMix = methylationMax+ minReverse
        p.patch(indexOrdered, methylationMix, alpha=0.5, line_width=2)
    p.y_range=Range1d(0,1.1)



#Figure definieren auf der alles abgebildet wird
p = figure(height=500, width=1000, name='plot')

#Restliche Widgets für Bedienung definieren

#DropDown Menü für Chromosomenauswahl, Start- und Endauswahl
def changes(attr, old, new):
    global startValue, endValue, chooseChromosom, checkboxRawAlgorithm, radioButtonSampleClass, finalGridplot, p
    #TODO überprüfen ob start und end auch wirklich Zahlen sind
    chromosom = str(chooseChromosom.value)
    start = int(startValue.value)
    end = int(endValue.value)
    if(finalGridplot.children[0].name == 'plot'):
        finalGridplot.children.remove(finalGridplot.children[0])
    else:
        print("This should not happen")
    p = figure(width=1000, height=500, name='plot')
    sample = True
    difference = False
    if(radioButtonSampleClass.active == 0):
        sample = True
        difference = False
    elif(radioButtonSampleClass.active == 1):
        sample = False
        difference = False
    # unterschied wurde ausgewählt
    else:
        sample = False
        difference = True
    # Sollen Rohdaten angezeigt werden
    if (0 in checkboxRawAlgorithm.active):
        displayDataRaw(chromosom, start, end, sample, difference)
    # Sollen AlgorithmenDaten angezeigt werden
    if(1 in checkboxRawAlgorithm.active):
        displayDataAlgorithm(chromosom, start, end, sample, difference)
    finalGridplot.children.insert(0, p)

#Auswahl zwischen Rohdaten, Klassen und Klassenunterschied
checkboxRawAlgorithm = CheckboxButtonGroup(
    labels=['Rohdaten', 'Algorithmus'], active=[0])
checkboxRawAlgorithm.on_change('active', changes)

radioButtonSampleClass = RadioButtonGroup(
        labels=['Sample', 'Klassen', 'Unterschied'], active=0)
radioButtonSampleClass.on_change('active', changes)

chromosoms = list(indexFile.keys())
chooseChromosom = Select(title='Chromosom', value="1", options=chromosoms)
chooseChromosom.on_change('value', changes)

startValue = TextInput(value='0', title="Start")
startValue.on_change('value', changes)

endValue = TextInput(value='100', title="Ende")
endValue.on_change('value', changes)

finalGridplot = layout(children = [p, [widgetbox(checkboxRawAlgorithm, radioButtonSampleClass, chooseChromosom), row(widgetbox(startValue, endValue))]], toolbar_location='left')
curdoc().add_root(finalGridplot)
