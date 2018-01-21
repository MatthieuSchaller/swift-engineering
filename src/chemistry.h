/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/
#ifndef SWIFT_CHEMISTRY_H
#define SWIFT_CHEMISTRY_H

/**
 * @file src/chemistry.h
 * @brief Branches between the different chemistry functions.
 */

/* Config parameters. */
#include "../config.h"
#include "chemistry_struct.h"

/* Import the right chemistry definition */
#if defined(CHEMISTRY_NONE)
#include "./chemistry/none/chemistry.h"
#elif defined(CHEMISTRY_GEAR)
#include "./chemistry/gear/chemistry.h"
#elif defined(CHEMISTRY_EAGLE)
#include "./chemistry/EAGLE/chemistry.h"
#else
#error "Invalid choice of chemistry function."
#endif

/* Common functions */
void chemistry_init(const struct swift_params* parameter_file,
                    const struct unit_system* us,
                    const struct phys_const* phys_const,
                    struct chemistry_data* data);

void chemistry_print(const struct chemistry_data* data);

#endif /* SWIFT_CHEMISTRY_H */