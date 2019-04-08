import sys
import numpy as np
from random import randint

import h5py

from bokeh.models.widgets import CheckboxButtonGroup, TextInput, Select
from bokeh.models.widgets import RadioButtonGroup, Div
from bokeh.io import show, curdoc
from bokeh.plotting import figure
from bokeh.layouts import layout, widgetbox, row, gridplot, column
from bokeh.models import LinearAxis, Range1d, HoverTool, ColumnDataSource

INTERVAL = 1000
SPLIT = " "
C = 1

def nearest_neighbour(array, first, second):
    """
    Findet bei einem sortierten Array das rechte erste Element und das linke
    zweite Element

    array (np.array): Array, indem die Werte gefunden werden sollen
    first (int) : Erste Wert, linke Schranke
    second (int) : Zweiter Wert, rechte Schranke

    """
    start = np.searchsorted(array, first, side="left")
    end = np.searchsorted(array, second, side="left")

    # Fall start und end außerhalb links
    if (first < array[0] and second < array[0]):
        cpg_start = 0
        cpg_end = 0
    # Fall start außerhalb end mittig
    elif (first < array[0] and 0 < end < len(array)):
        cpg_start = 0
        cpg_end = end
    # Fall start und end mittig
    elif (0 <= start < len(array) and 0 < end < len(array)):
        cpg_start = start
        cpg_end = end
    # Case start mittig end außerhalb
    elif (0 <= start < len(array) and end >= len(array)):
        cpg_start = start
        cpg_end = len(array)
    # Case start und end außerhalb
    else:
        cpg_start = len(array)-1
        cpg_end = len(array)-1
    return([cpg_start, cpg_end])

def sum_data(index, data, sample_name_list):
    """
    Addiert alle absoluten Daten nach der jeweiligen Klasse zusammen.

    Param:
        index (list): Liste von darzustellenden Indizes
        data (list): Liste der darzustellenden Daten
        sample_name_list (list): Namen der einzelnen Datenpunkte

    Returns:
        data_sumed (list): Zusammengefasste Daten nach Sample oder Klasse,
            data_sumed enthält die methylierten und unmethylierten Werte in
            der Form einer Liste von Dictionarys
            [{'sample_name', 'index', 'methylated', 'unmethylated'}]

    """

    # Bei gleichem Index wird methyliertes/ unmethyliertes zusammengerechnet
    data_sumed = []
    for x in range(len(data)):
        sample = sample_name_list[x]
        found = False
        found_index = -1
        methylation_methylated_sum = []
        methylation_unmethylated_sum = []
        sample_name = []
        # Sucht, ob entsprechendes Klasse schon vorhanden ist
        for i, dataset in enumerate(data_sumed):
            if sample in dataset['sample_name']:
                found = True
                found_index = i
        for i, dataset in enumerate(data[x]):
            # Finden der gleichen Indizes und addieren der Werte
            if found:
                previous_value = data_sumed[found_index]['methylated'][i]
                next_value = dataset[0] + dataset[2]
                new_value = previous_value + next_value
                data_sumed[found_index]['methylated'][i] = new_value
                previous_value = data_sumed[found_index]['unmethylated'][i]
                next_value = dataset[1] + dataset[3]
                new_value = previous_value + next_value
                data_sumed[found_index]['unmethylated'][i] = new_value
            # Sonst wurde ein neue Klasse gefunden und muss hinzugefuegt werden
            else:
                sample_name.append(sample)
                methylated = dataset[0] + dataset[2]
                unmethylated = dataset[1] + dataset[3]
                methylation_methylated_sum.append(methylated)
                methylation_unmethylated_sum.append(unmethylated)
        # Datensatz muss noch eingefuegt werden, wenn nicht gefunden
        if not found:
            data_values={
                'sample_name': sample_name,
                'index': index,
                'methylated':methylation_methylated_sum,
                'unmethylated':methylation_unmethylated_sum}
            data_sumed.append(data_values)
    return data_sumed


def sum_data_algorithm(index, data, class_name_list):
    """
    Soll alle Mittelwerte der Algorithmendaten bei gleichen Klassen berechnen.
	Aktuell gibt es nur Index, Daten und Klassennamen zurück.
	
    Param:
        index (list): Liste von darzustellenden Indizes
        data (list): Liste der darzustellenden Daten
        class_name_list (list): Namen der einzelnen Klassen

    Returns:
        index_sumed (list): Liste mit Index
        data_sumed (list): Liste mit Daten zusammenaddiert und
            gewichtet pro Klasse
        class_name_sumed (list): Liste der Namen
    """

    index_sumed = index
    data_sumed = data
    class_name_sumed = class_name_list
    return [index_sumed, data_sumed, class_name_sumed]

def make_data_sample(index, data, sample_name_list):
    """
        Fasst Rohdaten zu einem Datensatz zusammen.

        Param:
            index (list): Index der darzustellenden Punkte
            data (list) : Daten in einer mehrdimensionalen Liste
                            Es wird pro Sample die Daten an dem Index angegeben
            sample_name_list (list) : Namen pro Sample

        Returns:
            data_sumed (list) : Liste von Dictionarys,
                	jedes Element beschreibt ein Sample der Form
                    {'sample_name', 'index', 'methylated', 'unmethylated'}

    """

    data_sumed = []
    for x in range(len(data)):
        sample = sample_name_list[x]
        methylation_methylated_sum = []
        methylation_unmethylated_sum = []
        sample_name = []
        for y in range(len(data[x])):
            sample_name.append(sample)
            # Addieren der jeweils methylierten und unmethylierten Wert
            methylated = data[x][y][0] + data[x][y][2]
            unmethylated = data[x][y][1] + data[x][y][3]
            methylation_methylated_sum.append(methylated)
            methylation_unmethylated_sum.append(unmethylated)
        # Datensatz muss eingefuegt werden
        data_values={
            'sample_name': sample_name,
            'index': index,
            'methylated':methylation_methylated_sum,
            'unmethylated':methylation_unmethylated_sum}
        if data_values['index']:
            data_sumed.append(data_values)
    return data_sumed


def make_data_algorithm(index, data, names_sample_class):
    """
    Zusammenfassen der Daten des Algorithmus.

    Params:
        index (list) : Liste der darzustellendne Indizes
        data (list) : Daten mit einer Liste pro Sample, die eine
            Liste mit Methylierungslevel und Fehlerwert enthält
        names_sample_class (list) : Name der Samples oder Klassen

    Returns:
    data_sumed (list): Zusammengefasste Daten nach Sample oder Klasse,
        data_sumed enthält die methylierten und unmethylierten Werte in
        der Form einer Liste von Dictionarys
        [{'sample_name', 'index', 'methylation_level', 'error_value',
        'methylation_error_start', 'methylation_error_end'}]
    data_position (list): Enthält pro Sample jede Position der
        CpG-Dinukleotide

    """

    data_numpy = np.array(data)
    data_calculated = []
    for i, name in enumerate(names_sample_class):
        names = []
        index_list = []
        methylation_level_array = np.array([])
        error_list = []
        error_start_array = np.array([])
        error_end_array = np.array([])
        # Pro Sample Daten vorbereiten und konkatenieren
        for x, dataset in enumerate(data_numpy[i]):
            names.append(name)
            index_list.append(index[x])
            methylation_level_array = np.concatenate(
                (methylation_level_array, dataset[0]),
                axis=None)
            error_start_array = np.concatenate(
                (error_start_array, dataset[0] - dataset[1]),
                axis=None)
            for i, value in enumerate(error_start_array):
                if value < 0:
                    error_start_array[i] = 0
            error_end_array = np.concatenate(
                (error_end_array, dataset[0] + dataset[1]),
                axis=None)
            for i, value in enumerate(error_end_array):
                if value > 1:
                    error_end_array[i] = 1
            error_list.append(dataset[1])

        if names:
            data_calculated.append({
                'sample_name': names,
                'index': index_list,
                'methylation_level': methylation_level_array.tolist(),
                'error_value': error_list,
                'methylation_error_start': error_start_array.tolist(),
                'methylation_error_end': error_end_array.tolist()})
    return data_calculated


def calculate_methylation_level(data_sumed):
    """
    Berechnet die Methylierungslevel und Fehlerwerte von Datenpunkten.

    Params:
        data_sumed (list): Die Daten mit absoluten Werten von methylierten und
            unmethylierten Werten der einzelnen Samples. Sie haben die Form:
            [{'sample_name', 'index', 'methylated', 'unmethylated'}]

    Returns:
        data_list (list): Jeder Eintrag enthält Name, Index, Methyierungslevel,
            die absolute Anzahl an Methylierungen, den Fehlerwert und
            seinen Start- und Endwert
            [{'sample_name', 'index', 'methylation_level','methylated_absolut',
            'error_value', 'methylation_error_start', 'methylation_error_end'}]
        data_position (list): Jeder Eintrag enthält Name, Index und den Wert 0
            für die Darstellung der existierenden CpG-Dinukleotide

    """

    data_list = []
    data_position = []
    # Geht alle Eintraege der zusammengefassten Daten durch
    for data_sumed_sample in data_sumed:
        methylation_sum_list = []
        error = []
        start_error = []
        end_error = []
        methylation_level = []
        methylation_level_approx = []
        index = []
        index_position = []
        index_position_value = []
        methylation_position = []
        sample_name_data = []
        # Berechnung fuer die einzelnen Datenpunkte
        for y in range(len(data_sumed_sample['index'])):
            # Berechnung des Methylierungslevel
            methylation_sum = (data_sumed_sample['methylated'][y]
                + data_sumed_sample['unmethylated'][y])

            index_position.append(data_sumed_sample['index'][y])
            index_position_value.append(0)

            if methylation_sum != 0:
                sample_name_data.append(
                    data_sumed_sample['sample_name'][y])
                methylation_sum_list.append(methylation_sum)
                # Berechnung des Methylierungslevels für den Fehler
                meth_level_approx_calculated = (
                    data_sumed_sample['methylated'][y]
                    + 1 * C) / (methylation_sum + 2 * C)
                methylation_level_approx.append(meth_level_approx_calculated)

                index.append(data_sumed_sample['index'][y])
                # Berechnung des Erwartungswerts
                level = data_sumed_sample['methylated'][y] / methylation_sum
                methylation_level.append(level)

                # Berechnung der Standardabweichung
                sigma = (np.sqrt((meth_level_approx_calculated)
                    * (1-meth_level_approx_calculated))) / np.sqrt(
                    methylation_sum + 2 * C)
                error.append(sigma)
                if level-sigma < 0:
                    start_error.append(0)
                else:
                    start_error.append((level - sigma))
                if level+sigma > 1:
                    end_error.append(1)
                else:
                    end_error.append((level + sigma))
        data_values = {
            'sample_name': sample_name_data,
            'index': index,
            'methylation_level': methylation_level,
            'methylated_absolut' : methylation_sum_list,
            'error_value' : error,
            'methylation_error_start' : start_error,
            'methylation_error_end' : end_error}
        data_position_sample = {
            'sample_name': data_sumed_sample['sample_name'],
            'index_position': index_position,
            'index_position_value': index_position_value,
        }
        # Wenn leerer Eintrag, wird dieser nicht eingefuegt
        if data_values['sample_name']:
            data_list.append(data_values)
            data_position.append(data_position_sample)
    return data_list, data_position


def display_difference(methylation_level_all, data_calculated, difference, p):
    """
        Zeigt die Differenz an.

        Params:
            methylation_level_all (list): Alle Methylierungslevel pro Sample
                und Index
            data_calculated (list): Berechnete Daten entlang des Indizes
            difference (bool): True -> Differenz wird gezeigt,
                False -> Sie wird nicht gezeigt
            p (bokeh.plotting.figure.Figure): Die Figure, in der die Differenz
                angezeigt wird

    """

    # Wenn Differenz gefordert, muss diese noch ermittelt werden
    if difference:
        methylation_min = []
        methylation_max = []
        index_reverse = []
        min_reverse = []
        index_compared = data_calculated[0]['index']
        # Da es sich um unterschiedlich lange Daten handelt, muessen
        # diese auf einen gemeinsamen Index gebracht werden
        for index_sample in data_calculated:
            index_compared = list(set(index_compared).intersection(index_sample['index']))
        index_compared = sorted(index_compared)
        # Methylation_level_all auf dieses Intervall zuschneiden
        for i, methylation_level in enumerate(methylation_level_all):
            start_index = data_calculated[i]['index'].index(index_compared[0])
            end_index = data_calculated[i]['index'].index(index_compared[len(index_compared)-1])
            methylation_level_all[i] = methylation_level[start_index:end_index+1]

        # Finden von Minimum und Maximum, um Unterschied hervorzuheben
        methylation_min = list(np.amin(methylation_level_all, axis=0))
        methylation_max = list(np.amax(methylation_level_all, axis=0))
        # Daten werden in die richtige Reihenfolge fuer das Vieleck gebracht
        index_reverse = index_compared[:]
        index_reverse.reverse()
        index_ordered = index_compared + index_reverse
        min_reverse = methylation_min
        min_reverse.reverse()
        methylation_mix = methylation_max + min_reverse
        # Anzeigen des Unterschiedes
        p.patch(index_ordered, methylation_mix, alpha=0.5, line_width=2)


def display_data_raw(data_calculated, cpg_position, difference, line, error,
    absolute_bar, p, colors):
    """
    Zeigt die Rohdaten an.

    Params:
        data_calculated (list): Berechnete Daten mit dict{sampleName, index,
            Methylierungslevel, absolute Methylierungen, Fehler,
            Startwert des Fehlers und Endwert des Fehlers}
        cpg_position (list): Pro Sample die gefundenen CpG-Dinukleotide
            dict{'sample_name', 'index_position', 'index_value' }
        difference (bool): True -> Differenz soll angezeigt werden
            False -> Differenz soll nicht angezeigt werden
        line (bool): True -> Linie zwischen den Datenpunkten wird gezeichnet,
            False -> Keine Linie wird gezeichnet
        error (bool): True -> Fehlerwerte werden dargestellt,
            False -> Fehlerwerte werden nicht dargestellt
        absolute_bar (bool): True -> Absolute Werte werden als Balken dar-
            gestellt, False -> Absolute Werte werden nicht dargestellt
        p (bokeh.plotting.figure.Figure): Figure auf der die Daten angezeigt
            werden
        colors (list): Liste von Farben fuer die Darstellung im Diagramm

    """

    methylation_level_all = []
    # Sample oder Klassen durchgehen und Punkte, Balken und Linien zeichnen
    for i, data_calculated_dataset in enumerate(data_calculated):
        index_shift_points = [(value + i * 0.1)
            for value in data_calculated_dataset['index'][:]]
        index_shift_bars = [(value + i * 1)
            for value in data_calculated_dataset['index'][:]]
        data_calculated_dataset['index_shift_points'] = index_shift_points
        data_calculated_dataset['index_shift_bars'] = index_shift_bars
        cpg_position[i]['index_shift_points'] = [(value + i * 0.1)
            for value in cpg_position[i]['index_position'][:]]
        source_data = ColumnDataSource(data=data_calculated_dataset)
        source_position = ColumnDataSource(data=cpg_position[i])
        # Zeichnen der Punkte
        if error:
            error = p.vbar(
                x = 'index_shift_points',
                width = 0.01,
                bottom = 'methylation_error_start',
                top = 'methylation_error_end',
                color = colors[i],
                line_alpha = 0.1,
                fill_alpha = 0.1,
                source=source_data)
        points = p.x(
            'index_shift_points',
            'methylation_level',
            source = source_data,
            size = 10,
            color = colors[i],
            legend = data_calculated_dataset['sample_name'][0])
        position = p.cross(
            'index_shift_points',
            'index_position_value',
            source = source_position,
            size = 10,
            color = colors[i],
            line_alpha = 0.3,
            fill_alpha = 0.3,
            legend=cpg_position[i]['sample_name'][0])
        # Zeichnen der Linie
        if line:
            p.line(
                'index_shift_points',
                'methylation_level',
                source = source_data,
                color = colors[i],
                legend = data_calculated_dataset['sample_name'][0])
        # Zeichnen der Balken
        if absolute_bar:
            bars = p.vbar(
                x = 'index_shift_bars',
                width = 1,
                bottom = 0,
                top = 'methylated_absolut',
                color = colors[i],
                y_range_name = "absolute_data",
                fill_alpha = 0.10,
                line_alpha = 0.1,
                source = source_data)
            # Hovern ermoeglichen
            p.add_tools(HoverTool(
                renderers = [bars],
                tooltips = [
                    ("index", "@index"),
                    ("methylierung absolut", "@methylated_absolut"),
                    ("sample", "@sample_name")]))
        # Hovern ermoeglichen
        p.add_tools(HoverTool(
            renderers = [points],
            tooltips = [
                ("index","@index"),
                ("methylation_level","@methylation_level"),
                ("Fehler", "@error_value"),
                ("sample", "@sample_name")]))
        methylation_level_all.append(
            data_calculated_dataset['methylation_level'])

    display_difference(methylation_level_all, data_calculated, difference, p)
    if absolute_bar:
        # Berechnung des Maximum der absoluten Werte
        # Fuer eine bessere y-Achsen Darstellung
        max_absolut = 0
        for data_calculated_dataset in data_calculated:
            maximum_y_range = list(
                data_calculated_dataset['methylated_absolut'])
            if np.max(maximum_y_range) > max_absolut:
                max_absolut = np.max(maximum_y_range)
        # Hinzufuegen der Achsen
        p.extra_y_ranges = {"absolute_data":Range1d(
            start = 0, end = 2 * max_absolut)}
        p.add_layout(LinearAxis(y_range_name = "absolute_data"), 'right')

    p.y_range = Range1d(-0.01,1.1)

def display_data_algorithm(data_calculated, cpg_position, p, difference, error,
    colors, line):
    """
        Anzeigen der Daten des Algorithmus.

        Params:
            data_calculated (list): Pro Sample die zusammegefassten Werte
                in der Form dict{sampleName, index, Methylierungslevel,
                absolute Methylierungen, Fehler, Startwert des Fehlers und
                Endwert des Fehlers}
            cpg_position (list): Liste der CpG-Positionen
            p (bokeh.plotting.figure.Figure) : Figur auf der dargestellt wird
            difference (bool) : False -> Differenz wird nicht angezeigt,
                True -> Differnz wird angezeigt
            error (bool): True -> Fehler wird angezeigt,
                False -> Fehlerwerte werden nicht angezeigt
            colors (list) : Liste von Farben fuer die Samples oder Klassen
            line (bool) : False -> Keine Linie zwischen den Daten zeichnen,
                True -> Linie zwischen den Punkten zeichnen

    """

    methylation_level_all = []

    # Alle Sample oder Klassen durchgehen und Punkte zeichnen
    for i, dataset in enumerate(data_calculated):
        index_shift_points = [(value + i * 0.1)
            for value in dataset['index'][:]]
        dataset['index_shift_points'] = index_shift_points
        cpg_position[i]['index_shift_points'] = [(value + i * 0.1)
            for value in cpg_position[i]['index_position'][:]]
        source_data = ColumnDataSource(data=dataset)
        source_position = ColumnDataSource(data=cpg_position[i])

        if error:
            error = p.vbar(
                x = 'index_shift_points',
                width = 0.01,
                bottom = 'methylation_error_start',
                top = 'methylation_error_end',
                color = colors[i],
                line_alpha = 0.1,
                fill_alpha = 0.1,
                source = source_data)
        points = p.cross(
            'index_shift_points',
            'methylation_level',
            source = source_data,
            size = 10,
            color = colors[i],
            legend = dataset['sample_name'][0])
        position = p.cross(
            'index_shift_points',
            'index_position_value',
            source = source_position,
            size = 10,
            color = colors[i],
            line_alpha = 0.3,
            fill_alpha = 0.3)
        if line:
            p.line(
                'index_shift_points',
                'methylation_level',
                source = source_data,
                color = colors[i],
                legend = dataset['sample_name'][0])
        p.add_tools(HoverTool(
            renderers = [points], tooltips =[
                ("Index","@index"),
                ("Algorithmus Methylierungslevel","@methylation_level"),
                ("Fehler", "@error_value"),
                ("Sample", "@sample_name")]))
        methylation_level_all.append(dataset['methylation_level'])
    # Anzeigen der Differenz, wenn gewuenscht
    display_difference(methylation_level_all, data_calculated, difference, p)
    p.y_range = Range1d(-0.01,1.1)


def read_arguments():
    """
        Einlesen der Argumente nach dem Starten.

        Returns:
             [files_list, index_list, colors] : Eine Liste, die
                die Datei der Samples, die Indexdatei und zufällig
                generierte Farbwerte, enthält
				
		Raises:
			RuntimeError : Wenn nicht genau zwei Argumente übergeben werden

    """
    # Einlesen der Konfigurationsdatei und Abspeichern der einzelnen Dateien
    files_list = []
    colors = []
    if len(sys.argv) == 2:
        config_file = open(sys.argv[1], 'r')
        config_file_content = config_file.readlines()
        # Lesen der Zeilen: Startet bei der ersten Zeile, den Ueberschriften
        for line_config in config_file_content:
            # Sample, Class, Raw, Derived: 4 Spalten
            # Spalten bei dem gegebenen SPLIT Punkt
            line = line_config.split(SPLIT)
            # In files_list eintragen
            files_list.append(line)
    else:
        raise RuntimeError("Die korrekte Anzahl an Argumenten uebergeben "
            + "bokeh serve --show {Dateinamen} --args {Pfad zur Config Datei}")

    # Definieren der index_file
    index_file = h5py.File(files_list[1][2], 'r')
    files_list = files_list[2:]
    # Generieren von zufaelligen Farben, um diese spaeter anzuzeigen
    colors = [(randint(0,255), randint(0,255), randint(0,255))
        for i in range(len(files_list))]
    return [files_list, index_file, colors]


def prepare_data(chromosom, start, end, raw, sample, files_list, index_file,
    error_message):
    """
        Vorbereiten der Daten, bei gegebenem Chromosom, Start- und Endwert.

        Params:
            chromosom (str) : Darzustellendes Chromosom
            start (str) : Startwert der Daten
            end (str) : Endwert der Daten
            raw (bool) : True -> Rohdaten anzeigen,
                False -> Algorithmendaten anzeigen
            sample (bool) : True -> Sample anzeigen,
                False -> Klassen anzeigen
            files_list (list) : Liste von Samples mit allen Daten
            index_file (h5py._hl.files.File') : Indexdatei
            error_message (bokeh.models.widgets.markups.Div') :
                Widget fuer Fehlermeldungen in der Bedienung

        Returns:
            list bei Rohdatendarstellung : Pro Sample enthaelt diese ein
                Dictionary mit 'sample_name', 'index', 'methylated',
                'unmethylated'
            list bei Algorithmendarstellung :  Liste mit Index, Daten und
                Samplenamen

    """

    # Auslesen der aktuellen Werte
    try:
        start = int(start)
        end = int(end)
    except ValueError:
        error_message.text="Bitte eine Zahl eingeben"
        return
    # Ueberpruefe, ob das Intervall in den gegeben Grenzen liegt
    names_sample_class = []
    index = []
    if start >= 0 and end > 0 and start <= end and end - start <= INTERVAL:
        # Hole Daten von allen Samples
        # Start und End muessen definiert werden
        index = np.array(index_file.get(chromosom).get('cpg_positions'))
        # Berechnen der anzuzeigenden CpG-Dinukleotide
        cpg_intervall = nearest_neighbour(index, start, end)
        start_index = cpg_intervall[0]
        end_index = cpg_intervall[1]
        index = list(index[start_index: end_index])
        index_positions = []
        if not index:
            return
        data = []
        for i in files_list:
            index_positions.append(index)
            # Einlesen und Anhaengen der Methylierungsdaten
            try:
                if raw:
                    chromosom_file = h5py.File(i[2], 'r').get(chromosom)
                    dataset_all = chromosom_file.get('methylation')
                    dataset = dataset_all[start_index:end_index]
                else:
                    chromosom_file = h5py.File(i[3], 'r').get(chromosom)
                    firstDataset = list(chromosom_file.keys())[0]
                    dataset = chromosom_file.get(
                        firstDataset)[start_index:end_index]
                data.append(dataset)
            except OSError:
                error_message.text = "Die Datei konnte nicht geoeffnet werden."
                return
            except IndexError:
                text = "Bitte eine Zahl im gueltigen Bereich angeben."
                error_message.text= text
                return
            # Wenn Sample aktiv ist, pro Sample einlesen
            error_message.text = ""
            if sample:
                names_sample_class.append(i[0] + " " + i[1])
            else:
                names_sample_class.append(i[1])
    else:
        text = ("Start und Ende muessen groesser als 0 sein und " +
            "die Differenz zwischen diesen Zahlen muss kleiner als 1000 sein.")
        error_message.text = text
        return
    # Wenn Sample und Rohdaten angezeigt werden
    if sample and raw:
        return make_data_sample(index, data, names_sample_class)
    # Wenn Klasse und Rohdaten angezeigt werden
    elif not sample and raw:
        return sum_data(index, data, names_sample_class)
    # Wenn Algorithmendaten dargestellt werden sollen
    # Klassen und ALgorithmus, Daten muessen erst zusammengefasst werden
    elif not sample and not raw:
        return [index_positions,
            sum_data_algorithm(index, data, names_sample_class)]
    else:
        return [index_positions, [index, data, names_sample_class]]


def main():
    """
        Main-Methode, wird am Anfang aufgerufen.
    """
    # Daten aus den Argumenten auslesen
    dataFiles = read_arguments()
    files_list = dataFiles[0]
    index_file = dataFiles[1]
    colors = dataFiles[2]
    # Eingestellte Werte in der Anwendung
    chromosom = "1"
    start = "10400"
    end = "11000"
    raw_algorithm = [0]
    sample_class_diff = 0
    line = False
    error = True
    absolute_bar = True
    color_input_list = []

    # Figure definieren, auf der alles abgebildet wird
    p = figure(height=500, width=1000, name='plot')

    # Bei jedem Update der Widgets wird diese aufgerufen
    def start_routine():
        """
            Startet die passenden Methoden zum Veraendern der Grafik.
        """
        nonlocal chromosom, start, end, sample_class_diff, raw_algorithm, p

        # Entfernen der alten Figur
        if final_gridplot.children[0].name == 'plot':
            final_gridplot.children.remove(final_gridplot.children[0])
        # Anlegen einer neuen Figur, die die neuen Daten uebergeben bekommt
        p = figure(width=1000, height=500, name='plot')
        # setzen von sample und difference, je nachdem was gerade aktiv ist
        sample = True
        difference = False
        # sample ist aktiv, difference muss dann auf False gesetzt werden
        if sample_class_diff == 0:
            sample = True
            difference = False
        # Klasse ist aktiv
        elif sample_class_diff == 1:
            sample = False
            difference = False
        # Unterschied wurde ausgewaehlt, da dieser nur mit Klassen angezeigt
        # wird, muss sample ebenfalls auf False gesetzt werden
        else:
            sample = False
            difference = True
        # Bei Farbauswahl wird die Farbe dementsprechend angepasst
        for i, color_string in enumerate(color_input_list):
            color = color_string.value.split(",")
            try:
                color[0] = int(color[0])
                color[1] = int(color[1])
                color[2] = int(color[2])
            except ValueError:
                error_message.text = ("Bitte eine Zahl"
                    + " zwischen 0 und 255 eingeben.")
                return
            if (color[0]<0 or color[0]>255 or color[1]<0 or color[1]>255
                or color[2]<0 or color[2]>255):
                error_message.text = "Eine Zahl zwischen 0 und 255 eingeben."
                return
            colors[i] = (color[0], color[1], color[2])

        # Sollen Rohdaten angezeigt werden
        if 0 in raw_algorithm:
            sumed_data = prepare_data(chromosom, start, end, True, sample,
                files_list, index_file, error_message)
            if sumed_data is None:
                # Einfuegen der Figur
                final_gridplot.children.insert(0, p)
                return
            calculated_methylation = calculate_methylation_level(sumed_data)
            calculated_data = calculated_methylation[0]
            cpg_position = calculated_methylation[1]
            display_data_raw(calculated_data, cpg_position, difference,
                line, error, absolute_bar, p, colors)
        # Sollen Algorithmendaten angezeigt werden
        if 1 in raw_algorithm:
            pre_data = prepare_data(chromosom, start, end, False,
                sample, files_list, index_file, error_message)
            prepared_data = pre_data[1]
            cpg_position_samples = []
            # Dictionary fuer Index Positionen anlegen
            for index_sample in pre_data[0]:
                # Soll bei Null liegen
                index_value = [0 for i in index_sample]
                cpg_position = {'index_position': index_sample,
                    'index_position_value': index_value}
                cpg_position_samples.append(cpg_position)
            if prepared_data is None:
                # Einfuegen der Figur
                final_gridplot.children.insert(0, p)
                return
            calculated_data = make_data_algorithm(prepared_data[0],
                prepared_data[1], prepared_data[2])
            display_data_algorithm(calculated_data, cpg_position_samples,
                p, difference, error,
                colors, line)
        # Einfuegen der Figur
        final_gridplot.children.insert(0, p)

    # Restliche Widgets fuer Bedienung definieren
    def changes_raw_algorithm(attr, old, new):
        """
            Bei Veraenderungen bei der checkbox_raw_algorithm startet
            start_routine()

            Params:
                attr (str) : Gibt an, welches Attribut veraendert wurde
                old (bokeh.core.property.containers.PropertyValueList) :
                    Liste der alten Werte
                new (bokeh.core.property.containers.PropertyValueList) :
                    Liste der nun aktiven Werte

        """

        nonlocal raw_algorithm
        raw_algorithm = []
        for x in new:
            raw_algorithm.append(x)
        start_routine()

    def changes_sample_class_diff(attr, old, new):
        """
            Bei Veraenderungen bei der radio_button_sample_class startet
            start_routine()

            Params:
                attr (str) : Gibt an, welches Attribut veraendert wurde
                old (bokeh.core.property.containers.PropertyValueList) :
                    Liste der alten Werte
                new (bokeh.core.property.containers.PropertyValueList) :
                    Liste der nun aktiven Werte

        """

        nonlocal sample_class_diff
        sample_class_diff = new
        start_routine()

    def changes_chromosom(attr, old, new):
        """
            Bei Veraenderungen bei choose_chromosom startet start_routine()

            Params:
                attr (str) : Gibt an, welches Attribut veraendert wurde
                old (bokeh.core.property.containers.PropertyValueList) :
                    Liste der alten Werte
                new (bokeh.core.property.containers.PropertyValueList) :
                    Liste der nun aktiven Werte

        """

        nonlocal chromosom
        chromosom = new
        start_routine()

    def changes_start(attr, old, new):
        """
            Bei Veraenderungen bei start_value startet start_routine()

            Params:
                attr (str) : Gibt an, welches Attribut veraendert wurde
                old (bokeh.core.property.containers.PropertyValueList) :
                    Liste der alten Werte
                new (bokeh.core.property.containers.PropertyValueList) :
                    Liste der nun aktiven Werte

        """

        nonlocal start
        start = new
        start_routine()

    def changes_end(attr, old, new):
        """
            Bei Veraenderungen bei end_value startet start_routine()

            Params:
                attr (str) : Gibt an, welches Attribut veraendert wurde
                old (bokeh.core.property.containers.PropertyValueList) :
                    Liste der alten Werte
                new (bokeh.core.property.containers.PropertyValueList) :
                    Liste der nun aktiven Werte

        """

        nonlocal end
        end = new
        start_routine()

    def changes_line(attr, old, new):
        """
            Bei Veraenderungen von checkbox_button_line_circle startet
            start_routine()

            Params:
                attr (str) : Gibt an, welches Attribut veraendert wurde
                old (bokeh.core.property.containers.PropertyValueList) :
                    Liste der alten Werte
                new (bokeh.core.property.containers.PropertyValueList) :
                    Liste der nun aktiven Werte

        """

        nonlocal line
        if not new:
            line = False
        else:
            line = True
        start_routine()

    def changes_error(attr, old, new):
        """
            Bei Veraenderungen von checkbox_button_error startet
            start_routine()

            Params:
                attr (str) : Gibt an, welches Attribut veraendert wurde
                old (bokeh.core.property.containers.PropertyValueList) :
                    Liste der alten Werte
                new (bokeh.core.property.containers.PropertyValueList) :
                    Liste der nun aktiven Werte

        """
        nonlocal error
        if not new:
            error = False
        else:
            error = True
        start_routine()

    def changes_absolute_bar(attr, old, new):
        """
            Bei Veraenderungen von checkbox_button_error startet
            start_routine()

            Params:
                attr (str) : Gibt an, welches Attribut veraendert wurde
                old (bokeh.core.property.containers.PropertyValueList) :
                    Liste der alten Werte
                new (bokeh.core.property.containers.PropertyValueList) :
                    Liste der nun aktiven Werte

        """

        nonlocal absolute_bar
        if not new:
            absolute_bar = False
        else:
            absolute_bar = True
        start_routine()

    def changes_color(attr, old, new):
        """
            Startet die Routine, um die Farbwerte neu anzupassen
        """

        start_routine()

    # Auswahl zwischen Rohdaten, Klassen und Klassenunterschied
    checkbox_raw_algorithm = CheckboxButtonGroup(
        labels=['Rohdaten', 'Algorithmus'], active=[0])
    checkbox_raw_algorithm.on_change('active', changes_raw_algorithm)
    # Auswahl zwischen Sample und Klassen
    radio_button_sample_class = RadioButtonGroup(
            labels=['Sample', 'Klassen', 'Unterschied'], active=0)
    radio_button_sample_class.on_change('active', changes_sample_class_diff)
    # Auswahl fuer Liniendarstellung
    checkbox_button_line_circle = CheckboxButtonGroup(
            labels=["Linie"]
    )
    checkbox_button_line_circle.on_change('active', changes_line)
    # Auswahl fuer Fehlerdarstellung
    checkbox_button_error = CheckboxButtonGroup(
            labels=["Fehler"], active=[0]
    )
    checkbox_button_error.on_change('active', changes_error)
    # Auswahl fuer Darstellung der absoluten Daten
    checkbox_button_absolute = CheckboxButtonGroup(
            labels=["Absolute Werte"], active=[0]
    )
    checkbox_button_absolute.on_change('active', changes_absolute_bar)

    #DropDown Menue fuer Chromosomenauswahl, Start- und Endauswahl
    chromosoms = list(index_file.keys())
    choose_chromosom = Select(title='Chromosom', value="1", options=chromosoms)
    choose_chromosom.on_change('value', changes_chromosom)

    start_value = TextInput(value='10400', title="Start")
    start_value.on_change('value', changes_start)

    end_value = TextInput(value='11000', title="Ende")
    end_value.on_change('value', changes_end)

    # Fuer das Setzen der Farbwerte
    for i, color_value in enumerate(colors):
        color_input = TextInput(
            value= (str(color_value[0]) + "," + str(color_value[1]) + "," +
                str(color_value[2])), title="Farbwert " + str(i))
        color_input.on_change('value', changes_color)
        color_input_list.append(color_input)

    # Textfelder fuer Fehlermeldungen in der Anwendung
    error_message = Div(text="", width=300, height=10)

    final_gridplot = layout(children = [
        p,
        [column(widgetbox(
            checkbox_raw_algorithm,
            radio_button_sample_class,
            choose_chromosom)),
        widgetbox(checkbox_button_line_circle,
            checkbox_button_error, checkbox_button_absolute),
        row(
            widgetbox(start_value, end_value)),
        error_message],
        color_input_list],
        toolbar_location='left')
    curdoc().add_root(final_gridplot)

main()
