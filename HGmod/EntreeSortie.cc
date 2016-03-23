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

extern "C" {
#include <unistd.h>         // pour getopt
}

#include "HGmod.h"
#include "Messages.h"

using namespace std;
using namespace Utils_bg;

extern string program_name;
extern bool verbose;
extern bool vtk_format;
extern Messages msg;

void HGmod::print_usage (std::ostream& stream, int exit_code)
{
    stream << "\n *** HGmod - " << msg.getString("Un générateur de modèles hydrogéophysiques") << " ***\n\n";
    stream << msg.getString("Usage: ") << program_name << " [options] -p " << msg.getString("fichier.par") << '\n';
    stream << "  -h  " << msg.getString("Affiche ce message") << '\n'
    << "  -p  " << msg.getString("Spécification du fichier paramètres (obligatoire)") << '\n'
    << "  -v  " << msg.getString("Affiche les messages informatifs") << '\n'
    << "  -k  " << msg.getString("Sauvegarde des fichiers en format VTK")
    << std::endl;
    exit (exit_code);
}


string HGmod::traite_arg(int argc, char *argv[])
{
    string fichier_param;
    // -----------------------------------------------------------------
    //
    // Traitement des arguments
    //
    // -----------------------------------------------------------------
    
    int next_option;
    const char* const short_options = "hp:vk";
    bool pas_d_option = true;
    do {
        next_option = getopt (argc, argv, short_options);
        switch (next_option)
        {
                
            case  'h' : // -h or --help
                // User has requested usage information. Print it to standard
                //   output, and exit with exit code zero (normal termination).
                pas_d_option = false;
                print_usage (cout, 0);
                
            case  'p' :
                // This option takes an argument, the name of the input file.
                pas_d_option = false;
                fichier_param = optarg;
                break;
                
            case  'v' : // -v or --verbose
                pas_d_option = false;
                verbose = true;
                break;
                
            case 'k' :
#ifdef VTK
                pas_d_option = false;
                vtk_format = true;
#else
                std::cout << "\nWarning: VTK support not included during compilation.  Option -k ignored" << std::endl;
#endif
                break;
                
            case  '?' : // The user specified an invalid option.
                // Print usage information to standard error, and exit with exit
                //code one (indicating abnormal termination).
                print_usage (cerr, 1);
                
            case -1: // Done with options.
                break;
                
            default: // Something else: unexpected.
                abort ();
        }
    } while (next_option != -1);
    if ( pas_d_option || fichier_param == "" ) print_usage (cout, 0);
    
    return fichier_param;
}
