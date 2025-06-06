Pb: les intervalles proposés par CdfSolver n'étant pas optimaux, il n'est pas certains qu'il soit possible de trouver des potentiels pour les obtenir.
    En particulier, dans la version actuelle, rien ne garantie que les interfaces entre 2 agrégats soient bien positionnées avec les poids proposés.
    Une des solutions pourrait être de trouver les poids pour des interfaces au milieu des trous formés par les agrégats, et d'ajuster les poids
        Si on ajuste globalement, on pourrait chercher le min_w pour éviter les cellules vides et le max_w pour s'assurer qu'il y a au moins une cellule avec bords boule... mais rien ne garanti qu'on puisse toujours obtenir un interval non-vide.

Pour faire correctement l'initialisation, il faudrait résoudre un nouveau problème à chaque fois qu'on ajoute une nouvelle cellule à un agrégat
    Pour 1 cellule, on veut minimiser E( u ) = int( ( x - c )^2 * p( x ), x, I1, CINV( u + l / 2 ) )
        E'( u ) = ( I0 - c )^2 * p( I0 ) * CINV'( u + l / 2 ) -
                  ( I1 - c )^2 * p( I1 ) * CINV'( u + l / 2 )

        E'( u ) = I0 * ( I0 - 2 * c ) -
                  I1 * ( I1 - 2 * c )
        
        On cherche
                ( I0 - c )^2 = ( c - I1 )^2
        C'est à dire
                c = ( I0 + I1 ) / 2
        Ce qui signifie qu'on cherche u tel que c soit au milieu de I0 et de I1
        On pourrait construire CINV et faire une somme avec offset, qu'on inverse pour trouver la position optimale...
    Pour 2 cellules
        On a
            E( u ) = int( ( x - c0 )^2 * p( x ), x, CINV( u + 0 * l ), CINV( u + 1 * l ) )
                   + int( ( x - c1 )^2 * p( x ), x, CINV( u + 1 * l ), CINV( u + 2 * l ) )
            E'( u ) = ( I1( u ) - c0 )^2 - ( I0( u ) - c0 )^2
                    + ( I2( u ) - c1 )^2 - ( I1( u ) - c1 )^2
        Si on développe
            E'( u ) = I1( u )^2 - 2 * I1( u ) * c0 + c0^2
                    - I0( u )^2 + 2 * I0( u ) * c0 - c0^2
                    + I2( u )^2 - 2 * I2( u ) * c1 + c1^2
                    - I1( u )^2 + 2 * I1( u ) * c1 - c1^2
        Du coup:
            E'( u ) = 
                    + I2( u )^2 - 2 * I2( u ) * c1
                    - I0( u )^2 + 2 * I0( u ) * c0
                    + 2 * I1( u ) * ( c1 - c0 )
            E'( u ) = 
                    + I2( u ) * ( I2( u ) - 2 * c1 )
                    + I0( u ) * ( 2 * c0 - I0( u ) )
                    + 2 * I1( u ) * ( c1 - c0 )
        Si on dit I1( u ) = I0( u ) + l1( u ) et I2( u ) = I1( u ) + l2( u ) 
            E'( u ) = 2*I0*l1 + 2*I0*l2 - 2*c0*l1 - 2*c0*l2 - 2*d1*l2 + l1**2 + 2*l1*l2 + l2**2
            E'( u ) = I0*(2*l1 + 2*l2) + c0*(-2*l1 - 2*l2) - 2*d1*l2 + l1**2 + 2*l1*l2 + l2**2
            E'( u ) = 2*I0*(l1 + l2) - 2*c0*(l1 + l2) - 2*d1*l2 + (l1+l2)**2
            E'( u ) = ( 2*I0 - 2*c0 + l1 + l2 )*(l1 + l2) - 2*d1*l2
        On peut déjà dire que E'( u ) est croissante entre CDF( c0 ) et CDF( c1 )
        Si on met le bord gauche en c0, on a I0( u ) = c0. On peut proposer I1( u ) = c0 + o1 et I2( u ) = c0 + o2 avec o1 > 0 et o2 > o1
            E'( c0 ) = 
                     + ( c0 + o2 ) * ( ( c0 + o2 ) - 2 * c1 )
                     + c0 * ( 2 * c0 - c0 )
                     + 2 * ( c0 + o1 ) * ( c1 - c0 )
            E'( c0 ) = 
                     + c0^2 + c0 * o2 - 2 c0 * c1 + o2 * c0 + o2^2 - 2 * c1 * o2
                     + 2 * c0^2 - c0^2
                     + 2 * c0 * c1 - 2 * c0^2 + 2 * o1 * c1 - 2 * o1 * c0
            E'( c0 ) = 
                     + 2 * c0 * o2 + o2^2 - 2 * c1 * o2                    
                     + 2 * o1 * c1 - 2 * o1 * c0
            E'( c0 ) = 
                     + o2 * ( o2 - 2 * ( c1 - c0 ) )
                     + o1 * ( 2 * ( c1 - c0 ) )
            E'( c0 ) = o2^2 
                     + 2 * ( o1 - o2 ) * ( c1 - c0 )

        A priori, ce E'( u ) est une fonction croissante... On pourra toujours résoudre avec dichotomie "linéarisée"
        Cependant, pour le construire, on aura besoin de faire la somme de 3 fonctions pas définies dans le même espace
        Ou alors... on travaille sur des grilles régulière pour I( U ). Si les ratios sont constants, on pourrait avoir une grille d'un pas qui serait un ratio entier de `l`.
            E''( u ) = 
                    + I2'( u ) * ( I2( u ) - 2 * c1 ) + I2( u ) * I2'( u )
                    - I0'( u ) * ( I0( u ) - 2 * c0 ) - I0( u ) * I0'( u )
                    + 2 * I1'( u ) * ( c1 - c0 )
            E''( u ) = 
                    + 2 * I2'( u ) * I2( u ) - 2 * c1 * I2'( u )
                    - 2 * I0'( u ) * I0( u ) - 2 * c0 * I0'( u )
                    + 2 * I1'( u ) * ( c1 - c0 )
            E''( u ) / 2 = 
                    + I2'( u ) * ( I2( u ) - c1 )
                    + I0'( u ) * ( c0 - I0( u ) )
                    + I1'( u ) * ( c1 - c0 )

Rq: on pourrait ne pas se soucier de la borne supérieure pour les potentiels.
  Par exemple, on pourrait repartir sur l'idée de la convolution pour avoir une sensibilité sur les bords, 
  Si ratio est égale à 1, on utilise la CDF. On pourrait utiliser des inverses locaux dans les itérateurs pour gérer toutes les densités possibles
  Si le ratio est plus petit, ça serait avec une convolution

Partir de poids petits pour les faire grossir
  * on aurait des événements en grand nombre.


Pb de la sensibilité à la valeur absolue des potentiels
  * ça pourrait se regler simplement avec des convolutions
    * si ratio = 1, on fait tout avec CDF
    * sinon, on fait une convolution pour avoir toujours une sensibilité

Calcul des dérivées par rapport à la largeur de convolution
  * on cherche Aire calculée = aire prescrite pour chaque cellule
  * Rq: on aura quand même des "changements de régime" quand on passera des balles aux interfaces
    Est-ce qu'on pourrait couper des cellules avec des bords flous ?
    par exemple, on dirait que la densité d'une cellule, c'est la multiplication de fonctions step "douces"

Prop: on part d'une ICDF avec ratio = 1 et on fait une frontière floue pour les boules.

Dérivées:
  * M x = V 
    M[ i, j ] = -rho( x_j ) * dist_diracs si i != j et ( rho( x_j ) + rho( x_i ) ) * dist_diracs
    V[ i ] = integral( rho, x_i, x_j )
    x_j = (d_(j-1)+d_j)/2 + (w_j-w_(j-1)) / (d_j-d_(j-1))/2
  * M' x + M x' = V'
    M'[ i, j ] = -( x_j' * drho/dx( x_j ) + drho/dc( x_j ) ) * dist_diracs
      x_j' = (w_j'-w_(j-1)') / (d_j-d_(j-1))/2
    V[ i ] = x_i' * rho( x_i ) - ... + integral( drho/dc, x_i, x_j )

Ccl: c'est qd même un peu relou de trouver les dérivées...
Rk: de toutes façons, ça ne résoud pas le pb de l'insensibilité à la valeur absolue des potentiels s'ils sont trop hauts... Sauf si on met de la convolution au dela du support

Autre rk: dans le fond, on a un pb de contrôle.
  Lorsqu'on avance trop vite dans les potentiels, on se retrouve avec une direction du newton qui change
  On pourrait mettre une pénalisation non linéaire sur
  * les cellules vides à cause des différences de potentiel (par exemple en se basant sur une taille min, calculée avec le max de la densité)
  * les cellules vides à cause des potentiels négatifs
  * les cellules en dehors du domaine ou dans 
  A priori, la position des interfaces étant localement une fonction linéaire du potentiel, on peut trouver une limite.

Rq: on peut envisager de dilater des intervalles pour décrire les densités. Ça pourrait marcher pour les interfaces

On aimerait un coût qui empêche de se retrouver avec des dérivées nulles
  Le critère le plus naturel serait de se baser sur log( aire / aire_prescrite ), mais ça exclut les aires négatives
  Rq: le log va faire encore plus exploser les cellules sur les densités faibles.
    -> on pourrait aussi pénaliser la taille des cellules indépendamment de la densité, par exemple en se basant sur la max de la densité, qui définit une taille min de cellule.
    Le point, c'est que dans la solution, il n'est pas possible d'avoir une aire_densité beaucoup plus petite que aire_reelle : on sait que l'aire réelle ne peut pas être plus petite que la masse prescrite / le max de la densité.
      On aimerait par exemple un blend qui tend vers l'aire réelle lorsque ça devient plus petit que le min.

      Autre prop: on pondère par la densité moyenne
   
      densite_moyenne * log( aire / aire_prescrite )
    
  Rq: la non linéarité devrait donner un chemin non-linéaire. On peut évaluer les dérivées

  Rq: 

Considérant que A = V + M * w

On peut chercher à résoudre une EDP d( V + M * w )/dt = log( ( V + M * w ) ./ P )
  Si on injecte w = w1 * t + w2 * t^2 + w3 * t^3 + ... on obtient
    M * ( w1 + 2 * w2 * t + 3 * w3 * t^2 + ... ) = log( ( V + M * ( w1 * t + w2 * t^2 + w3 * t^3 + ... ) ) ./ P )
  en t = 0, on obtient
    M * w1 = log( V ./ P )
  En dérivant 1 fois
    M * ( 2 * w2 + 6 * w3 * t + ... ) = M * ( w1 + 2 * w2 * t + 3 * w3 * t^2 + ... ) * P ./ ( V + M * ( w1 * t + w2 * t^2 + w3 * t^3 + ... ) )
  en t = 0, on obtient
    M * ( 2 * w2 ) = M * w1 * P ./ V
  En dérivant 1 fois supplémentaire
    M * ( 6 * w3 + ... ) = M * ( 2 * w2 + 6 * w3 * t + ... ) * P ./ ( V + M * ( w1 * t + w2 * t^2 + w3 * t^3 + ... ) ) -
                           M * ( w1 + 2 * w2 * t + 3 * w3 * t^2 + ... ) * P ./ ( V + M * ( w1 * t + w2 * t^2 + w3 * t^3 + ... ) )^2
  en t = 0, on obtient
    6 * M * w2 = M * P * ( 2 * w2 / V - w1 / V^2 ) * P
    24 * M * w3 = M * ( 6 * w3 /  - 2 * w2 / V - w1 / V ) * P / V

-> on pourrait faire une descente locale de gradient à pas fixe, et on extrapole.

En 1D, on a qd même bien envie d'utiliser les CDF
  Qu'est-ce que ça pourrait donner de faire une transition vers le partiel ?

  M Xs = V
  M X' = V' - M' X avec X = Xp + a * Xs
  M'' X + 2 * M' X' + M X'' = V'' avec X = Xp + a * Xs + a^2/2 * Xd

Si on fait la convolution seule, on aura toujours le problème des coupes de boules.

Pb de la mort: les intersections avec les boules introduisant des non-linéarités, on se retrouve avec des problèmes de relaxation.
  Pire encore: on peut se retrouver avec des changement de régime...
  Si on oublie ces changements de régime, on pourrait.

Proposition: on cherche une convolution qui assure que ce n'est pas la densité qui va empêcher de converger.
  Si on oublie le Newton, on peut dire que ce qu'on veut éviter, c'est qu'une cellule devienne vide à cause d'un contraste trop fort.
  

--------------------------------------------------------------------------------------------------------------
Objectif une bonne semaine en musique, pratique collective.
Age des participants -> ensemble vocal + percussion. Supervision autonomie. Ce qui est propoposé c'est le cadre.
Feu de camp, boum.
Animations, manuel, créatif, dessins jeux de sociétés.
Adultes séparés pour groupes, chambres.
Chambres de 3, plus resseré. Elle sera pas seule dans sa tranche d'âge. + de 8 à 12. + de filles 60%.

Arrivée le dimanche ? Activités. 1 tiers la veille. Petit confort. Pas de nécessité. Possible seule dans la chambre.

Liens avec les autres styles ?

Ambiance de travail.

  45 enfants, principalement entre 7 et 13. 2 de 13 ans, 4 de 12 ans. Pas de 14 ans.
  Activités en groupes en fonction des ages. Ensemble dans le centre.

Impro 
Bernard

Séjour linguistique sur campus en angleterre. Archer. 35 de 10 à 16 ans. 1425e.
Dossier avec Sonia.
--------------------------------------------------------------------------------------------------------------

La proposition, c'est de controler la descente des densités. En gros on va controler l'amplitude des déplacements.
  Pour les bords extérieurs, il suffira

--------------------------------------------------------------------------------------------------------------


Dimanche soir possible.
Libéré le jeudi.
  23km. marche beaucoup. presque de la ballade.
  Prise en charge des bagages. Optionnel.
  Déjeuner pas inclus. Si option sandwich, ou sinon, endroit supérette, restaurant.

--------------------------------------------------------------------------------------------------------------
Puy en Velay -> 4h27, train à 19h00 gare de Lyon (1 changement)
Lannion -> 3h40
Colmar -> 2h23
Brest -> 3h51

Déjà fait chambon, pavin

sac a dos 30L, au plus light, se changer le soir. Le pantalon 5 jours, tongs pour le soir.

Prop: on cherche pour chaque ensemble connecté l'erreur en fonction 

--------------------------------------------------------------------------------------------------------------
Prop: plutôt que fixer les poids, on fixe la position des interfaces. Par exemple, on donne la position en cours pour bi et la position avec un offset
  Ça permet d'obtenir deux listes de positions, pour lesquelles l'erreur est nulle pour les cellules ii, et avec une position pour ib.
  Ensuite, on peut paramétrer en w0 : pour w0 fixé, la masse impose une frontière pour bi et donc une frontière pour ib puis un w  
  Si on impose w0, trouver la frontière demande à résoudre un pb non linéaire. D'où l'idée d'imposer plutôt la frontière.
  Du coup, on ne connait pas w0, mais on peut calculer w1 - w0

Que stocke-t-on ?
  * l'offset précédant des masses 

On peut aussi paraméter par l'ajout de masse sur la première interface.
  L'ajout de masse sur la dernière, c'est - l'ajout de masse sur la première

Rq de la mort: on pourrait simplement calculer la masse manquante entre les couped de boules.
  Pour un w_bi donne, on va donc pouvoir trouver de façon immédiate
    * le w_ib utile pour trouver la masse manquante (en utilisant l'integrale bord à bord)
    * la position 

Pour un w_0 donné,
  * on veut masse_prescrite = masse_actuelle + ( sqrt( w0 ) - sqrt( w0_o ) ) * density->value( b0 ) + écart_de_masse_à_droite
  * on est donc capables de trouver l'écart_de_masse_à_droite
  * et donc, on peut déterminer w1... jusqu'à wn

Prop: on cherche la masse de la dernière cellule comme polynôme de r0, le rayon de la première cellule.
  à chaque étape, on pourra donner le rayon fonction de r0

Prop: on contruit les rayon (racines des poids) fonction du rayon de BI.

Pb: les jonctions "non anticipées" de cellules semblent provoquer une augmentation importante du nombre d'itérations

Prop: on fait une initialisation plus poussée.
  L'idéal, ça serait de calculer les poids exacts à chaque ajout de cellule.
    Pb: les polynômes qui donnent les poids optimaux dépendent de la position dans la densité.
    Il faudrait donc réévaluer tous les poids à chaque changement de position.
  Rq: on pourrait faire une initilisation pour une densité constante et converger vers la densité finale.
  Autre prop: pour l'initialisation, on reste sur l'idée de supposer que les densités ne bougent pas, et en plus, on calcule en direct les poids
  Prop pour démarrer : on calcule les poids initiaux par dichotomie (on pourra optimiser après)

Prop 1: au moment de donner les poids, on optimize le poids initial pour trouver le poids final... -> on va tomber sur les mêmes pb de convergence
Prop 2: on fait le parcourt pour trouver les w dans les 2 sens, et on fait un blend
  On est certains de trouver des w positifs... mais est-ce qu'on est sûrs d'éviter les cellules vides ?
Prop 3: on fait une répartition égale des bords de cellule
  -> ce n'est pas une bonne façon de pré-initialiser, en particulier lorsqu'il y a des contrastes dans la densité
Prop 4: plutôt que de chercher beg_x et end_x, on cherche des w...
Prop 5: on fait un passage en partant de w0. Quand on arrive à la fin, on regarde de combien il y a eu plantage au niveau de la masse et on refait un passage avec la masse corrigée.
  Prop: on calcul la masse entre x0 et end_x

Prop: on fait un vecteur de polynômes qui donnent les interfaces en fonction du poids de la première cellule.
  L'objectif, c'est d'optimiser la poids de la première cellule à chaque ajout dans le pack.
  Pb: on peut donner la suite de poids uniquement pour une densité connue. Lorsqu'on décale les cellules, les densité peuvent changer radicalement et la représentation polynomiale devient dégueulasse.
  
TODO:
* dichotomie si Newton échoue
* vérifier les masses relatives
* gérer les diracs placés de façon identique
* vérifier si le poids final de l'init peut poser problème
* max_mass_ratio_error_target
* readme

Pb: la dérivée de l'erreur de positionnement des agrégats n'est pas une fonction croissante. Du coup, on ne peut pas vraiment utiliser les valeurs extrêmes pour sortir d'un Newton...
  
En gros, l'initialisation doit être parfaite pour pouvoir démarrer les newtons.
  Le problème, c'est qu'on a plein d'itérations lorsqu'on fonctionne en mode "conservatif"

Rq: si la densité est constante, il est trivial d'initialiser les poids, y compris en partiel.
  L'idée de regulariser la densité est séduisantes, mais ça veut dire qu'il faut calculer les dérivées de la solution.
  Quel critère prendre pour 

M x= V
M x' + M' x = V' (avec x = 0)
M x'' + 2 * M' x' + M'' x = V'' (avec x = 0)
M x''' + 3 * M' x'' + 3 * M'' x' + M''' x = V''' (avec x = 0)

M x''' + 4 * M' x''' + 6 * M'' x'' + 4 * M''' x' + M'''' x = V'''' (avec x = 0)

La proposition, c'est de tester une projection en regardant si le newton est capable de converger

( ( 1 - t ) I + t ∇ ) X = ( 1 - t ) D

min ∑ ( 2 * u_n - u_(n-1) - u_(n+1) )^2 / 2
-> ∑ ( 2 * u_n - u_(n-1) - u_(n+1) )

Le problème, c'est que ni Neuman ni Dirichlet ne fonctionnent, à moins de trouver comment les fixer 
  On pourrait par exemple faire une fonction de la valeur de départ, avec l'idée de minimiser une fonction quadratique du genre somme( ( u_n - u_(n-1) )^2 )...
  C'est quand même bien moins simple qu'un filtre. Ça ne va pas non plus bien se traduire en 2D ou +.

Rq: dans le fond, on voudrait une diffusion sur un domaine plus large pour ne garder qu'une partie.

Autre prop: on fait un blend entre la fonction et sa moyenne, mais on va plus ou moins vite pour y aller en fonction de la distance à un "évènement"... bof

Autre prop: on résoud avec Dirichlet de moins en moins local... Rq: c'est peu coûteux en 1D mais après, ça revient à faire n^1/d convolutions. 

Pb de la mort : les dérivées sont très discontinues. Prop: on essaye avec extrapolation

d3 - 3 * d2 + 2 * d1 - d0

M x = V
M x' + M' x = V' -> 

M x'' + 2 * M' x' + M'' x = V''



d = X * ( 1 / sum( x[...] ) )
d' = X' / sys_div - X * sum( x' ) / sum( x )^2

d'' = X'' / sys_div
    - 2 * X' * sum( x' ) / sum( x )^2
    + 2 * X * sum( x' )^2 / sum( x )^3
    - X * sum( x'' ) / sum( x )^2

Est-ce qu'on pourrait faire une initialisation avec epsilon très grand ? Ça serait une solution pour simplifier les calculs au démarrage“


mean( sr.time_in_solve ):  0.00504353
mean( sr.time_in_solve ):  0.00265965
Sliced transport ICP: 0.000611775 seconds per iteration
