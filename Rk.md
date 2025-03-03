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
                    - I0( u ) * ( I0( u ) - 2 * c0 )
                    + 2 * I1( u ) * ( c1 - c0 )
        A priori, ce E'( u ) est une fonction croissante... On pourra toujours résoudre avec dichotomie "linéarisée"
        Cependant, pour le construire, on aura besoin de faire la somme de 3 fonctions pas définies dans le même espace
        Ou alors... on définie les fonctions sur des grilles, et ça simplifie pas mal de choses
