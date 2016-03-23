/*
 *  Granulo.cc
 *  HGmod
 *
 *  Created by Bernard Giroux on 08-01-15.
 *  Copyright 2008 Bernard Giroux. All rights reserved.
 *
 */

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
#include "Granulo.h"

using namespace std;

void Granulo::calculFractionVolumique()
{
    vector<double> volume(fractionVolumique);
    double fractionMassique;
    double somme = 0.0;
    double sommeVol = 0.0;
    for (size_t n=0; n<fractionVolumique.size(); ++n) {
        fractionMassique = 1.0 - passant[n] - somme;
        somme += fractionMassique;
        volume[n] = fractionMassique / masseVolumique[n];
        sommeVol += volume[n];
    }
    for (size_t n=0; n<volume.size(); ++n) {
        fractionVolumique[n] = volume[n] / sommeVol;
    }
}