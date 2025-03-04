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
        Si on dit I1( u ) = I0( u ) + o1( u ) et I2( u ) = I0( u ) + o2( u ) 
            E'( u ) = 
                    + I2( u ) * ( I2( u ) - 2 * c1 )
                    + I0( u ) * ( 2 * c0 - I0( u ) )
                    + 2 * I1( u ) * ( c1 - c0 )
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

Prop intermédiaire: si 