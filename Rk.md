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