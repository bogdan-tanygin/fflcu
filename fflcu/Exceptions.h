/*
  Copyright (C) 2012 Alexander (Polyakov) Peletskyi

  This file is part of FFLCU.

  FFLCU is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  FFLCU is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

/*
 * Exceptions.h
 *
 *  Created on: Jan 31, 2012
 *      Author: alexander
 */

#ifndef EXCEPTIONS_H_
#define EXCEPTIONS_H_
#include <string>
	struct WrongDeviceNumberException {
		WrongDeviceNumberException (int n) {
			deviceNumber = n;
		}
		unsigned int deviceNumber;
	};

	struct TooDenseSystemException{
		TooDenseSystemException(unsigned int lastPart){
			lastParticleNumber = lastPart;
		}
		unsigned int lastParticleNumber;
	};

	struct DeviceMemoryException{
		DeviceMemoryException(std::string message){
			this->message = message;
		}
		std::string message;
	};

	struct HostMemoryAllocationException{
		HostMemoryAllocationException(std::string message){
			this->message = message;
		}
		std::string message;
	};

	struct DeviceMemoryAllocationException: public DeviceMemoryException{
		DeviceMemoryAllocationException(std::string message):DeviceMemoryException(message){}
	};

	struct DeviceMemoryCopyException: public DeviceMemoryException{
		DeviceMemoryCopyException(std::string message):DeviceMemoryException(message){}
	};

	struct DeviceMemCpyToSymbolException: public DeviceMemoryException{
		DeviceMemCpyToSymbolException(std::string message):DeviceMemoryException(message){}
	};

	struct DeviceMemsetException: public DeviceMemoryException{
		DeviceMemsetException(std::string message):DeviceMemoryException(message){}
	};

#endif /* EXCEPTIONS_H_ */
