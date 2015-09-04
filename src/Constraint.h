/*=========================================================================

* @file
* @author  Lin Ma <majcjc@gmail.com>
* @version 1.0
*
* @section LICENSE
*
* This program is free software; you can redistribute it and/or
* modify it under the terms of the GNU General Public License as
* published by the Free Software Foundation; either version 2 of
* the License, or (at your option) any later version.
*
* This program is distributed in the hope that it will be useful, but
* WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
* General Public License for more details at
* http://www.gnu.org/copyleft/gpl.html
*
* @section DESCRIPTION
*
* Constraint for fast mass spring.

=========================================================================*/

#ifndef Constraint_H
#define Constraint_H

class Constraint
{
public:
  Constraint();
  virtual ~Constraint();

  virtual void init() = 0;
  virtual void buildMatrix() = 0;
  virtual void buildRightHand() = 0;
  virtual void updateRightHand() = 0;
  virtual void updateProjection() = 0;

private:
  Constraint(const Constraint&); // not implemented
  void operator = (const Constraint&); // not implemented
};

#endif