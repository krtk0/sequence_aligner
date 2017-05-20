from Substitution_Matrix import pam120, pam250, blosum62
from Aligners import Global_Aligner, Local_Aligner
from appJar import gui

""" GUI for Aligner. ONGOING """

app = gui("Aligner", "500x300")

app.addEntry("Sequence #1")
app.addEntry("Sequence #2")

app.setEntryDefault("Sequence #1", "SEQUENCEONEHERE")
app.setEntryUpperCase("Sequence #1")
app.addButton("Submit 1", )

app.setEntryDefault("Sequence #2", "SEQUENCETWOHERE")
app.setEntryUpperCase("Sequence #2")

app.go()