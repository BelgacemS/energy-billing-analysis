# Plan du projet : erreurs de facturation et prévision de consommation des sites mobiles

## Objectif

Analyser les factures d'énergie mensuelles de 92 sites d'accès radio (35 mois, janvier 2019 à novembre 2021) pour répondre à deux questions :

- Q1 : quelles factures paraissent erronées, pourquoi, avec quel écart en kWh ? S'agit-il d'une erreur ponctuelle ou d'un changement durable du site ?
- Q2 : quelle consommation prévoir en décembre 2021 pour chaque site, avec quelle marge d'erreur ? Sur quels sites cette prévision est-elle peu fiable ?

Le livrable est un rapport d'analyse de style entreprise, sans code, lisible par un responsable non technique.

### Critère de réussite global

Depuis un clone propre du dépôt :

1. `Rscript -e 'renv::restore()'` puis `Rscript run_all.R` régénèrent toutes les tables, figures et chiffres du rapport dans `outputs/`, sans intervention manuelle.
2. `Rscript -e 'testthat::test_dir("tests/testthat")'` passe sans échec.
3. `quarto render report/rapport.qmd` produit un PDF qui contient :
   - les statistiques descriptives et leurs enseignements ;
   - une liste priorisée de factures suspectes, chacune avec son motif, son type, son écart estimé en kWh et le sens de cet écart ;
   - le taux de détection et le taux de fausses alertes de la méthode, mesurés par injection d'anomalies ;
   - une prévision de décembre 2021 pour les 92 sites, avec des intervalles à 80 % et 95 % dont la couverture a été mesurée sur des mois jamais utilisés pour les choix ;
   - la liste des sites peu fiables, avec deux scénarios chiffrés pour ceux dont la fin de série est ambiguë.
4. Aucun chiffre du rapport n'est recopié à la main : tous sont lus dans `outputs/results/` par du code masqué dans le rendu.

## Contraintes retenues

- Données : un seul fichier, 3 220 factures (colonnes `NAME`, `BEGIN`, `END`, `KWH_PERIOD`), aucune donnée externe.
- Pas de vérité terrain sur les erreurs : la détection est non supervisée et doit être argumentée.
- La prévision s'appuie d'abord sur les méthodes du cours : méthodes simples (naïve, saisonnière naïve), lissage exponentiel, ARIMA/SARIMA, validation croisée temporelle, diagnostic des résidus.
- Rapport sans code ni sortie brute de logiciel : chiffres, tableaux et graphiques commentés, méthodes expliquées en langage clair.
- Ton du rapport : celui d'un data scientist qui présente son étude à un service achats ou énergie. Rédaction sobre et directe, conclusions chiffrées, limites dites franchement.
- Travail en binôme, pas de soutenance. Niveau visé : complet.
- Le fichier de données est versionné tel quel dans `data/` (dépôt public partagé par le binôme) et n'est jamais modifié : tout nettoyage se fait dans le code. Les documents du sujet restent hors du dépôt.

## Choix techniques

- R avec `forecast` : ce sont les fonctions du cours (`tsCV`, `monthdays`, `ses`, `ets`, `auto.arima`, `snaive`, `accuracy`). C'est aussi l'implémentation de référence d'auto.arima (Hyndman & Khandakar, 2008).
- `renv` pour l'environnement : versions des paquets figées dans `renv.lock`.
- Unité d'analyse, le log des kWh par jour calendaire du mois : retire l'effet de longueur du mois (février vaut 0,91 de la médiane annuelle du site en brut, 0,98 par jour, comme janvier et mars) et rend les écarts comparables entre sites dont la facture médiane va de 544 à 433 000 kWh par mois. On divise par les jours du mois et non par les jours facturés, pour ne pas fabriquer de faux pic sur une période partielle : une période irrégulière est un signal de Q1 traité par les contrôles.
- Effet du mois estimé en commun sur les sites (médianes) : la saison réelle est faible (par jour, −1,5 à −1,7 % en hiver et +3 % en août) et l'effet propre à chaque site se reproduit mal d'une année sur l'autre (ρ de Spearman entre −0,3 et +0,4, médiane ≈ 0,1) : estimé sur deux ans, il serait surtout du bruit, d'où la mise en commun (Withycombe, 1989). Les 2 gros sites gardent leur propre saisonnalité.
- Détection par panel robuste : niveau du site découpé en paliers par segmentation PELT (Killick et al., 2012, paquet `changepoint`, variante non paramétrique dans `changepoint.np`), plus effet commun du mois, résidus normalisés par une échelle robuste. La segmentation sépare ruptures et pics, là où une médiane glissante prend les bords d'une rupture pour des pics. Or 74 % des sites ordinaires (67 sur 90) ont un changement de niveau de plus de 10 %.
- Filtre de Hampel par site (médiane glissante, MAD) gardé comme méthode de référence : il chiffre le gain du panel sur le même banc d'essai.
- Typologie à règles explicites : chaque facture signalée s'explique en une phrase à un responsable non technique.
- Seuils, pénalité de segmentation et coût choisis par injection d'anomalies simulées, avec une règle de choix écrite avant le premier lancement : seul moyen de mesurer détection et fausses alertes sans vérité terrain, sans régler les seuils sur ce qu'on veut trouver.
- Analyse d'intervention de Chen et Liu (1993) écartée : fragile sur 35 points avec plusieurs ruptures. Le rapport la cite comme approche classique.
- Backtest à origine glissante : validation sur les prévisions de juillet 2020 à mai 2021 (11 origines, dont décembre 2020), test sur juin à novembre 2021 (6 origines) lu une seule fois. Les choix se font sur la validation, la précision annoncée vient du test.
- Règle de sélection de la méthode de prévision écrite avant de lancer : MASE moyenne de validation la plus basse, la méthode la plus simple si l'écart est de moins de 2 %, choix par classe et jamais site par site (11 points de validation par site sont trop peu).
- MASE (Hyndman & Koehler, 2006) normalisée par la naïve par jour : comparable entre sites, sans division par zéro.
- Intervalles empiriques (quantiles des erreurs relatives de validation, mis en commun par classe de fiabilité) comparés aux intervalles gaussiens des modèles : on garde la construction dont la couverture de test est la plus proche du nominal.
- Rapport en Quarto (`.qmd` rendu en PDF par Typst, code masqué) : format standard des rapports de data science, chiffres et figures produits par le code au moment du rendu.

## Organisation du dépôt

```
PLAN.md  README.md  run_all.R  renv.lock  .gitignore
data/            # données brutes (versionné, jamais modifié)
R/               # fonctions : config, chargement, contrôles, détection, prévision
scripts/         # une étape d'analyse par script, appelés par run_all.R
tests/testthat/  # tests des fonctions critiques
outputs/         # figures, tables, results/ (non versionné, régénéré)
report/          # rapport.qmd
```

Chemins et paramètres (seuils, pénalité, origines du backtest, règle de sélection, graines aléatoires) centralisés dans `R/config.R`, chemin des données modifiable par la variable d'environnement `DATA_PATH`.

## Étapes

### Mise en place

- [x] **Étape 1 : initialiser le projet**
  - Objectif : dépôt git, environnement reproductible, squelette qui tourne.
  - Fichiers : `.gitignore` (`.DS_Store`, `outputs/`, `.env`, secrets, `.Rhistory`, `.RData`, `renv/library/`, PDF de sujet, `.quarto/`, PDF rendu), `.gitattributes` (fins de ligne LF, le binôme travaille sous Windows), `renv.lock`, `.Rprofile`, `R/config.R`, `run_all.R`, `report/rapport.qmd` (squelette, figures en PNG à 300 dpi : sous macOS, le SVG de R dépend de Cairo et donc de XQuartz), `report/.Rprofile` (active renv quand Quarto est lancé depuis `report/`, par exemple depuis un IDE ; depuis la racine, c'est le `.Rprofile` racine qui est lu), `README.md`. Le CSV est déplacé dans `data/`. Remote GitHub configuré.
  - Outils : Quarto (`brew install --cask quarto`) ; paquets R `forecast`, `changepoint`, `changepoint.np`, `ggplot2`, `dplyr`, `tidyr`, `testthat`, `knitr`, `rmarkdown`.
  - Vérification : `Rscript run_all.R` affiche le nombre de lignes lues (3 220) ; `quarto render report/rapport.qmd` produit un PDF d'une page qui affiche ce nombre ; `git status` ne montre que des fichiers à versionner (le CSV de `data/` y figure, ni PDF du sujet ni sorties).
  - Hors périmètre : tout traitement des données.

- [x] **Étape 2 : charger et normaliser les données**
  - Objectif : table longue propre (site, suffixe, mois, début, fin, jours facturés, jours du mois, kWh, kWh/jour) et matrice 92 × 35 du log des kWh par jour calendaire.
  - Fichiers : `R/load.R`, `tests/testthat/test-load.R`.
  - Vérification : tests sur valeurs connues (0007BR janvier 2019 : 2 500 kWh, 31 jours ; 0628BR septembre 2021 : 16 jours facturés, 30 jours du mois) ; 92 sites × 35 mois, aucune valeur manquante après conversion du séparateur de milliers.
  - Hors périmètre : toute décision sur les anomalies.

### Statistiques descriptives

- [x] **Étape 3 : statistiques descriptives et enseignements**
  - Objectif : décrire niveaux, dispersion, tendance et saisonnalité. Établir aussi les constats qui justifient les choix : poids des 2 gros sites (environ 80 % des kWh), effet calendrier contre saison réelle, profil du mois commun ou non (globalement, par suffixe), reproductibilité de l'effet propre à chaque site d'une année sur l'autre, part des sites avec un changement de niveau, précision de la naïve brute contre la naïve par jour.
  - Contenu : distribution des niveaux médians (2 gros sites à part), profil mensuel brut contre par jour, dispersion des effets du mois entre sites, croissance 2019 à 2021 par site, carte de chaleur site × mois de l'écart au niveau du site, comparaison des trois suffixes.
  - Fichiers : `scripts/01_descriptive.R`, `R/features.R`, sorties dans `outputs/figures/`, `outputs/tables/`, `outputs/results/`.
  - Vérification : les ordres de grandeur du sujet sont retrouvés (médiane d'un site d'environ 2 260 kWh/mois, 2 sites au-dessus de 300 000 kWh) ; chaque figure a une phrase de conclusion notée dans le script.
  - Hors périmètre : détection et prévision.

- [ ] **Étape 4 : pipeline minimal de bout en bout**
  - Objectif : une première réponse complète aux deux questions avec les méthodes les plus simples, pour que le projet ait un résultat à chaque étape suivante.
  - Contenu : contrôles déterministes (périodes hors mois calendaire, trous et chevauchements, valeurs nulles ou négatives, doublons, valeurs identiques deux mois de suite comparées à la fréquence attendue par hasard) ; détection par filtre de Hampel ; prévision naïve par jour de décembre 2021 (× 31) ; squelette du rapport qui affiche les deux listes.
  - Fichiers : `R/checks.R`, `R/hampel.R`, `scripts/02_baseline.R`, `tests/testthat/test-checks.R`, `tests/testthat/test-hampel.R`, `report/rapport.qmd`.
  - Vérification : tests sur de petites tables construites à la main ; sur les données, 2 périodes irrégulières, 2 trous de couverture, 15 valeurs répétées sur 12 sites et leur fréquence attendue par hasard ; 92 prévisions sans valeur manquante ; le rapport squelette se génère.
  - Hors périmètre : panel, segmentation, modèles.

### Q1 : erreurs de facturation

- [ ] **Étape 5 : valeur attendue par panel et ruptures de niveau**
  - Objectif : pour chaque facture, une valeur attendue = palier du site + effet commun du mois. Les paliers viennent d'une segmentation PELT du log des kWh par jour corrigé de l'effet du mois, avec un segment minimal d'un mois pour qu'un pic forme son propre segment. Résidus normalisés par une échelle robuste du bruit (MAD des différences premières).
  - Fichiers : `R/panel.R`, `R/segment.R`, `scripts/03_panel.R`, `tests/testthat/test-segment.R`.
  - Vérification : sur séries synthétiques, une rupture franche est trouvée au bon mois, un pic isolé ne crée pas de rupture durable, une série constante bruitée ne produit aucune rupture ; sur les données, la rupture de 1240BR en mai 2021 est trouvée. Point à trancher dans cette étape : le coût gaussien avec `minseglen = 1` sépare-t-il vraiment pics et ruptures, en le comparant à `changepoint.np::cpt.np` sur les mêmes séries synthétiques ? Si aucun des deux n'y arrive, appliquer le repli décrit dans les risques avant de passer à l'étape 6.
  - Hors périmètre : choix définitif de la pénalité et du coût (étape 7).

- [ ] **Étape 6 : typologie des factures suspectes et écarts en kWh**
  - Objectif : classer chaque facture signalée parmi : période irrégulière, écart ponctuel (isolé, retour au palier), paire compensée (sous-facturation puis rattrapage, somme des écarts proche de 0), écart temporaire (quelques mois puis retour), changement durable (nouveau palier maintenu), valeur implausible, indéterminé (rupture trop proche de la fin). Pour chacune : écart en kWh par rapport à la valeur attendue, avec sens (sur ou sous-facturation) et fourchette, cumul depuis la rupture pour un changement durable. Classement par écart en kWh dans chaque niveau de confiance, 2 gros sites dans un tableau séparé.
  - Fichiers : `R/classify.R`, `scripts/04_suspects.R`, `tests/testthat/test-classify.R`, `outputs/tables/suspects.csv`.
  - Vérification : tests sur des séries synthétiques de chaque type ; les cas cités dans l'énoncé (0089KO mai à août 2019, 1067KO juillet 2019, 1088KO, 1978KO et 1986KO en 2020, 0254BB) et 1240BR apparaissent avec un type cohérent, ou leur absence est expliquée.
  - Hors périmètre : calibration.

- [ ] **Étape 7 : validation par injection d'anomalies et calibration**
  - Objectif : mesurer le comportement du détecteur sans vérité terrain. Construire des panels simulés à partir des données (paliers et effet du mois estimés, résidus rééchantillonnés), y injecter des anomalies de types et d'amplitudes connus (pics de 1 % à 100 %, paire compensée, écart temporaire, rupture de niveau), mesurer taux de détection, taux de fausses alertes et exactitude de la typologie. Comparer le panel au filtre de Hampel sur ce banc. Choisir pénalité, coût (gaussien ou non paramétrique) et seuils par la règle écrite dans `R/config.R` avant le premier lancement. Analyse de sensibilité : nombre d'alertes réelles quand les seuils varient.
  - Fichiers : `R/simulate.R`, `scripts/05_calibrate.R`, `tests/testthat/test-simulate.R`.
  - Vérification : courbes détection contre amplitude par type, pour les deux méthodes ; 20 réplications à graines fixées, résultats identiques à deux exécutions ; liste des suspects régénérée avec les paramètres retenus.
  - Hors périmètre : nouvelles méthodes de détection.

### Q2 : prévision du mois suivant

- [ ] **Étape 8 : banc d'essai sans fuite et méthodes de référence**
  - Objectif : backtest à origine glissante, prévision à un mois, validation de juillet 2020 à mai 2021, test de juin à novembre 2021. Tout prétraitement (effet du mois, segmentation, neutralisation des écarts ponctuels, classes, quantiles) est recalculé à chaque origine avec les seules données disponibles. Méthodes de référence : naïve brute, naïve par jour, saisonnière naïve. Métriques : MASE, erreurs relatives (médiane, moyenne, 90e centile), erreur totale en kWh, avec et sans les 2 gros sites. Règle de sélection écrite dans `R/config.R`.
  - Fichiers : `R/backtest.R`, `scripts/06_baselines.R`, `tests/testthat/test-backtest.R`.
  - Vérification : test de non-fuite (modifier les valeurs postérieures à une origine ne change pas la prévision faite à cette origine) ; erreurs des références sur la validation.
  - Hors périmètre : modèles du cours ; aucune lecture des origines de test.

- [ ] **Étape 9 : modèles du cours site par site**
  - Objectif : lissage exponentiel simple, ETS et auto.arima sur le log des kWh par jour de chaque site, en validation. Repérer où auto.arima fonctionne mal (ordres dégénérés, différenciation causée par une rupture, résidus autocorrélés au test de Ljung-Box) et expliquer pourquoi.
  - Fichiers : `R/models.R`, `scripts/07_models.R`.
  - Vérification : temps d'exécution mesuré d'abord sur 5 sites ; tableau des erreurs de validation par méthode ; liste des sites où auto.arima échoue, avec le motif.
  - Hors périmètre : classes et variantes robustes.

- [ ] **Étape 10 : classes, variantes robustes et sélection**
  - Objectif : classes de comportement par classification hiérarchique (Ward) sur des caractéristiques calculables à chaque origine (niveau, pente, volatilité, ruptures), 2 gros sites à part, nombre de classes choisi par dendrogramme, silhouette et interprétabilité, comparaison avec le suffixe. Candidats supplémentaires : modèle par classe sur la série corrigée de l'effet du mois mutualisé, niveau local robuste (lissage simple après neutralisation des écarts ponctuels et repartant du dernier palier), combinaison simple des meilleurs candidats (Bates & Granger, 1969). Application de la règle de sélection par classe.
  - Fichiers : `R/classes.R`, `R/robust.R`, `scripts/08_select.R`.
  - Vérification : tableau des classes (effectif, profil, exemples) ; tableau méthodes × classes sur la validation ; méthode retenue par classe et application de la règle tracée dans `outputs/results/`.
  - Hors périmètre : lecture du test.

- [ ] **Étape 11 : évaluation finale, marges d'erreur et prévisions de décembre 2021**
  - Objectif : lire une seule fois les origines de test avec la stratégie retenue. Intervalles à 80 % et 95 % : quantiles empiriques des erreurs relatives de validation mis en commun par classe de fiabilité, comparés aux intervalles gaussiens, construction retenue selon la couverture de test. Indicateur de fiabilité (méthode pas meilleure que la naïve, rupture ou anomalie en fin de série, intervalle large) validé : les sites signalés doivent avoir une erreur de test plus forte. Deux scénarios chiffrés (écart ponctuel ou changement durable) pour les sites dont la dernière facture est ambiguë. Prévisions finales réajustées sur les 35 mois.
  - Fichiers : `R/intervals.R`, `scripts/09_forecast.R`, `tests/testthat/test-intervals.R`, `outputs/tables/forecast_dec2021.csv`.
  - Vérification : couverture de test des deux constructions, écart au nominal discuté ; erreur de test des sites signalés contre non signalés ; 92 lignes, bornes ordonnées autour de la prévision.
  - Hors périmètre : horizon supérieur à un mois.

### Rapport

- [ ] **Étape 12 : sorties finales et reproductibilité**
  - Objectif : figures, tableaux et résultats définitifs dans `outputs/`, `run_all.R` qui enchaîne tous les scripts, README avec les commandes exactes.
  - Fichiers : `run_all.R`, `scripts/*.R`, `README.md`.
  - Vérification : après suppression de `outputs/`, `Rscript run_all.R` régénère tout ; temps total noté dans le README ; tests verts.
  - Hors périmètre : rédaction.

- [ ] **Étape 13 : rédaction du rapport**
  - Objectif : rapport de style entreprise : synthèse d'une page, données et hypothèses, statistiques descriptives, Q1 (méthode, justification, validation par injection, liste priorisée), Q2 (méthode, protocole, prévisions et marges, sites peu fiables et scénarios), limites et recommandations (obtenir les tarifs pour chiffrer en euros), annexes (liste complète des suspects, table des 92 prévisions).
  - Fichiers : `report/rapport.qmd` (document principal), qui assemble avec `{{< include >}}` une partie par fichier : `report/_summary.qmd`, `report/_data.qmd`, `report/_descriptive.qmd`, `report/_q1.qmd`, `report/_q2.qmd`, `report/_limits.qmd`, `report/_appendix.qmd`. Chaque membre du binôme peut rédiger sa partie en parallèle sans conflit git.
  - Vérification : `quarto render report/rapport.qmd` produit le PDF ; aucun code ni sortie brute visible ; chaque figure commentée ; chaque affirmation appuyée par un chiffre calculé ou une référence ; hypothèses de travail écrites explicitement.
  - Hors périmètre : nouvelle analyse.

## Risques et détection précoce

- Fuite de données dans le backtest (prétraitement calculé sur la série entière) : test de non-fuite dès l'étape 8.
- Sélection optimiste des seuils ou des méthodes : règles de choix écrites avant lancement, test lu une seule fois à l'étape 11.
- Absence de vérité terrain pour Q1 : l'injection mesure la détectabilité des motifs injectés, pas le taux réel d'erreurs. Le banc descend jusqu'à des amplitudes proches du bruit (environ 1 %) pour ne pas être trop facile. Le rapport parle de factures à vérifier, pas d'erreurs avérées.
- Confusion entre rupture et pic : test dédié à l'étape 5, comparaison au filtre de Hampel à l'étape 7.
- Coût de segmentation non robuste aux pics. Dans `changepoint`, PELT (`cpt.mean`) n'accepte que les coûts `"Normal"` (gaussien) et `"CUSUM"`. Aucun coût l1 ou robuste n'est disponible, alors qu'il en existe en Python (paquet `ruptures`). Deux parades sont prévues : un segment minimal d'un mois (`minseglen = 1`), pour qu'un pic forme son propre segment au lieu de déplacer un palier ; une comparaison à l'étape 7 avec la variante non paramétrique `changepoint.np::cpt.np` (PELT avec un coût fondé sur la distribution empirique). Ces deux API ont été vérifiées dans la documentation du CRAN. Leur comportement réel sur nos séries reste inconnu jusqu'à l'étape 5, où il se constate. Repli si les deux options confondent pics et ruptures : segmentation robuste aux valeurs aberrantes (Fearnhead & Rigaill, 2019), dont l'implémentation R reste à trouver, ou retrait des pics détectés avant la segmentation.
- Effet du mois mal estimé ou non partagé : visible à l'étape 3 (dispersion entre sites et entre suffixes). Repli : effet par groupe, ou aucun si négligeable en kWh par jour.
- Fins de série ambiguës (1240BR en tête) : type « indéterminé » en Q1, deux scénarios en Q2.
- Ruptures simultanées sur plusieurs sites (regroupements en juillet et octobre 2019, mai et juin 2020 ; juillet 2019 est aussi le premier mois testable et recueille des changements antérieurs) : plus probablement des changements du réseau que des erreurs de facturation. À l'étape 6, une rupture partagée par plusieurs sites le même mois n'est pas classée comme erreur sans autre indice.
- Classement par kWh écrasé par les 2 gros sites (1 % de 0698BR correspond à environ 4 300 kWh, presque deux mois d'un site médian) : tableau séparé dès l'étape 6.
- Temps de calcul du backtest (92 sites × 17 origines × plusieurs modèles) : mesure sur 5 sites avant le passage complet.
- Installation des paquets R et rendu Quarto : vérifiés dès l'étape 1 avec `renv::restore()` et un rapport squelette.

## Annexe A : éléments retenus du cadrage

- Questions posées : (1) y a-t-il des erreurs de facturation ? (2) comment prévoir la consommation du mois suivant pour chaque site ? Rapport attendu : statistiques descriptives et enseignements, justification des méthodes et conclusions pour les deux questions.
- Hypothèses de travail à écrire dans le rapport : `KWH_PERIOD` est l'énergie facturée entre `BEGIN` et `END`, en kWh, même unité pour toutes les lignes ; une période normale couvre un mois calendaire complet et un écart de dates est une anomalie de facturation ; un écart fort et ponctuel est plus probablement une erreur, un écart durable plus probablement un changement d'équipement ; le suffixe du nom est une variable de regroupement possible dont le sens est inconnu.
- Pistes évoquées pour la prévision : modèles par classes de comportement, ou auto-ARIMA par site avec repérage des sites où il fonctionne mal. Les deux sont évaluées comme candidates (étapes 9 et 10). Aucune indication pour la facturation.
- Premières observations à confirmer : périodes irrégulières (0628BR septembre 2021, 0089KO février 2020), sous-facturation suivie d'un rattrapage (0089KO 2019), ruptures de niveau simultanées (1088KO, 1978KO, 1986KO en 2020 ; 1067KO juillet 2019), perte de saisonnalité (0254BB à partir de 2020), 15 valeurs identiques consécutives sur 12 sites, forte tendance des 2 gros sites.

## Annexe B : constats vérifiés au cadrage

- 3 220 factures, 92 sites × 35 mois, aucune valeur manquante ni nulle ; exactement 2 périodes hors mois calendaire, qui laissent 2 trous de couverture.
- Les 2 gros sites représentent environ 81 % des kWh facturés.
- Sur les années complètes 2019 et 2020, février vaut 0,91 de la médiane annuelle du site en brut et 0,98 par jour, au niveau de janvier et mars : l'effet février est un effet calendrier. Par jour, l'écart saisonnier est de −1,5 à −1,7 % en hiver, +2,2 % en juillet et +3,2 % en août.
- La naïve à un mois a une erreur relative médiane de 3,7 % en brut et de 1,2 % par jour (de février 2019 à mai 2021, mois de test exclus) ; la naïve par jour est meilleure sur 70 % des factures.
- L'effet propre de chaque site se reproduit mal d'une année sur l'autre : ρ de Spearman entre −0,32 et +0,40 sur les 23 paires de mois et les deux centrages (par année civile ou glissant), médiane 0,06 ou 0,12 selon le centrage ; pour décembre 2019 et 2020, −0,09 ou 0,24. La persistance est faible et seulement probable : 17 corrélations positives sur 23 avec le centrage glissant (Wilcoxon p = 0,01), 13 sur 23 avec le centrage annuel (p = 0,10).
- 74 % des sites ordinaires (67 sur 90) ont un changement de niveau de plus de 10 % entre deux fenêtres de 6 mois ; les mois de rupture se regroupent (juillet et octobre 2019, mai et juin 2020 ; juillet 2019 est aussi le premier mois testable, ce qui gonfle son effectif). Les différences entre suffixes ne sont pas significatives.
- 1240BR passe d'environ 3 000 à environ 140 kWh par mois à partir de mai 2021 (baisse de 95 %), cas absent des observations de l'énoncé.
- Les 15 valeurs identiques deux mois de suite ne sont pas plus fréquentes qu'attendu par hasard : un écart nul apparaît 15 fois, un écart donné de 1 à 10 kWh environ 12,5 fois en moyenne. Elles ne suffisent pas à conclure à une facturation estimée.

## Références

- Bates, J. M. & Granger, C. W. J. (1969). The combination of forecasts. Operational Research Quarterly.
- Chen, C. & Liu, L.-M. (1993). Joint estimation of model parameters and outlier effects in time series. JASA.
- Fearnhead, P. & Rigaill, G. (2019). Changepoint detection in the presence of outliers. JASA.
- Hyndman, R. J. & Athanasopoulos, G. Forecasting: Principles and Practice (2e éd.), ouvrage support du cours.
- Hyndman, R. J. & Khandakar, Y. (2008). Automatic time series forecasting: the forecast package for R. Journal of Statistical Software.
- Hyndman, R. J. & Koehler, A. B. (2006). Another look at measures of forecast accuracy. IJF.
- Killick, R., Fearnhead, P. & Eckley, I. A. (2012). Optimal detection of changepoints with a linear computational cost. JASA.
- Makridakis, S. & Hibon, M. (2000). The M3-Competition: results, conclusions and implications. IJF.
- Makridakis, S., Spiliotis, E. & Assimakopoulos, V. (2020). The M4 Competition: 100,000 time series and 61 forecasting methods. IJF.
- Withycombe, R. (1989). Forecasting with combined seasonal indices. IJF.
