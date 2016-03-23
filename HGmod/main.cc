// Copyright (c) 2016 Bernard Giroux. All rights reserved.
/*
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 */
#include <iostream>
#include <fstream>

#include "HGmod.h"
#include "Grille.h"
#include "MilieuPoreux.h"
#include "Messages.h"

using namespace std;
using namespace HGmod;

string program_name;
bool verbose = false;
bool vtk_format = false;
Messages msg("fr");

int main (int argc, char *argv[]) {
	
    program_name = argv[0];
    
    string fichier_param = traite_arg(argc, argv);
    paramEntree<double> pe;
    get_param( fichier_param, pe );
	
    // Génère la grille
    Grille grl(pe.pGrille, pe.pSortie); 

    Regional* cov = new Regional();
    Generateur<Regional>* gen = 0;
    
    if ( pe.simulerPorosite ) {

        cov->remplire(pe.pPhys.pCov);
	
        switch ( pe.pGen.typeGen ) {
            case Utils_bg::TBM :
                gen = new TBM<Regional>(pe.pGrille, cov, pe.pGen.seed,
                                        pe.pGen.L);
                break;
			
            case Utils_bg::FFTMA :
                gen = new FFTMA<Regional>(pe.pGrille, cov, pe.pGen.seed);
                break;
			
            default:
                cerr << msg.getString("err1") << endl;
                abort();
        }
    }
    else {
        grl.lirePorosite(pe.fichierPorosite);
    }
    
    MilieuPoreux<double,Generateur,Regional> milieu( pe.pPhys, gen );
    
	if ( verbose ) {
		cout << "\n"<< msg.getString("Milieu Poreux")<< "\n";
		milieu.affiche_param( cout );
	}
	
    if ( !vtk_format )
        grl.sauveGrille(pe.pSortie.fichierOut);
	
	if ( pe.pPhys.sData.conditionne ) {
		grl.remplireSalinite(milieu);
        if (pe.pSortie.sauveSalinite)
            grl.sauveSalinite(pe.pSortie.fichierOut);
    }
    
    if ( pe.pGen.pCond.conditionne && pe.simulerPorosite )
        grl.prepareCondPorosite(milieu, pe.pGen.pCond);
    
	for (int n=1; n<=pe.nSimulations; n++) {
		if ( verbose )
            cout << "\n\n***" << msg.getString(" Simulation no ") << n << " ***\n";
		if (n!=1)
            gen->reinitialise();
		grl.remplire(milieu, pe.pGen.pCond, pe.f, pe.frequences);
        if ( vtk_format )
            grl.sauveVTK(pe.pSortie.fichierOut, n, pe.frequences);
        else
            grl.sauveData(pe.pSortie.fichierOut, n, pe.frequences);
	}
	
	delete gen;
	delete cov;
	
	cleanup( pe );
	
	return 0;
}
