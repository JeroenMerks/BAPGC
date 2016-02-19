Wat is dit?
--------------------------

Dit is een pipeline voor het ontdekken van genen buiten een opgegeven pathway
waarvan sterke aanwijzingen bestaan dat deze door dezelfde transcriptiefactoren
 worden gereguleerd. De gevonden genen worden op verschillende manieren verder onderzocht.


Benodigdheden
--------------------------

* Ubuntu Linux >= 14.04
* 4Gb of meer werkgeheugen
* Een internetverbinding


Installatie
--------------------------

Er is een Bash-script aanwezig die automatisch benodigde softwarepakketten
downloadt en installeert van de Ubuntu en Python pip repositories. Dit script wordt
 automatisch uitgevoerd na het starten van de pipeline.
Pas op! Het kan voorkomen dat je het wachtwoord van je gebruikersaccount (of su) dient op te geven!

Hier volgt een lijst van afhankelijkheden, mocht er toch een pakket blijken te missen:
* python 2.7.x
* build-essential (Ubuntu C compiler)
* Biopython >= 1.66 (Python tools for computational molecular biology)
* clustalw = 2.1
* pandas >= 17.1 (data analysis / manipulation library for Python)
  - NumPy >= 1.7.1 (package for scientific computing with Python)
  - setuputils (library designed to facilitate packaging Python projects)
* scipy
* mygene >= 2.3
* pylab
* ete3
* matplotlib >= 1.5.1


Gebruiken
--------------------------
Gebruik:
python(2) start_pipeline.py [PATHWAY_CODE]

PATHWAY_CODE, de KEGG pathway code van de metabole route

Voorbeeld:
python start_pipeline.py 'hsa04150'

Resultaten zijn te vinden in de /data/results/ folder.


Credits
--------------------------
De makers van MOODS en Melissa van Wieringen, Jeroen Merks, Koen van Diemen, Rick de Graaf en Christian van der Niet.