// Copyright 2016 Khaled Y. M. (https://github.com/kmehrunes)
/*
This file is part of GLex Genome Compression.

GLex is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License version 3
as published by the Free Software Foundation.

GLex is distributed in the hope that it will be useful, but 
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with GLex.  If not, see https://www.gnu.org/licenses/gpl.txt.
*/

using System;
using System.IO;

namespace GLex
{
	/// <summary>
	/// A reader for reading FASTQ files.
	/// </summary>
	public class FASTQReader
	{
		StreamReader reader;

		public FASTQReader (string file)
		{
			reader = new StreamReader (file);
		}

		/// <summary>
		/// Reads  a chunk of 4 lines where the first and the third are discarded regardless
		/// of the content. The second and forth lines represent the sequence and quality
		/// scores respectively.
		/// </summary>
		/// <returns>
		/// An array of length two of strings, where the first string is the sequence
		/// and the second one is the quality scores of each nucleotides
		/// </returns>
		public string[] ReadChunk ()
		{
			if (reader.EndOfStream)
				return null;
			
			try {
				string[] data = new string[2];

				// skip sequence ID
				reader.ReadLine ();
				// read the sequence
				data [0] = reader.ReadLine ();
				// skip separator 
				reader.ReadLine ();
				// read quality scores
				data [1] = reader.ReadLine ();

				if (data[1] == null)
					throw new Exception ("Failed to read quality scores of sequence " + data[0]);

				return data;
			}
			catch {
				return null;
			}
		}
	}
}

