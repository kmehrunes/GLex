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
	/// A writer class for GLex files, currently only supports DNA and only raw A,C,G,T strings.
	/// </summary>
	public class GLexWriter
	{
		StreamWriter writer;

		// information to be written to the header
		byte segmentLength = 0;
		ulong sequenceLength = 0;
		byte bits;

		// buffer-related variables
		char[] buffer;
		int buffered = 0;

		// used just to keep track of the segments while debugging
		int writtenSegments = 0;

		// Sequence length is used to create the header of the file (to know the length of the last segment)
		/// <summary>
		/// Initializes a new instance of the <see cref="GLex.GLexWriter"/> class.
		/// </summary>
		/// <param name="file">The path to the file to which the compressed genome will be written to</param>
		/// <param name="mode">The encoding mode of the file, which determines the number of bits used for encoding.</param>
		/// <param name="totalLength">
		/// The length of the full genome sequence. Used in the header, and to know the length of the last segment
		/// </param>
		public GLexWriter (string file, EncodingMode mode, ulong totalLength)
		{
			writer = new StreamWriter (file);
			segmentLength = (byte) GLexEncoding.MaxSequenceLength (mode);
			buffer = new char[segmentLength];
			sequenceLength = totalLength;

			switch (mode) {
			case EncodingMode.UNSIGNED_16_BITS:
				bits = 16;
				break;

			case EncodingMode.UNSIGNED_32_BITS:
				bits = 32;
				break;

			case EncodingMode.UNSIGNED_64_BITS:
				bits = 64;
				break;
			}

			WriteHeader ();
		}

		/// <summary>
		/// Writes the header to the file. Must be called before
		/// writing anything else to the file.
		/// </summary>
		private void WriteHeader ()
		{
			byte[] bytes;

			writer.BaseStream.WriteByte (bits);

			bytes = BitConverter.GetBytes (sequenceLength);
			writer.BaseStream.Write (bytes, 0, bytes.Length);

			writer.BaseStream.WriteByte (segmentLength);
		}

		/// <summary>
		/// Writes the given sequence to the file. The sequence does not have to be
		/// of a certain length, the writer handles buffering and encoding the
		/// segments.
		/// </summary>
		/// <param name="sequence">The sequence of a pattern (a segment)</param>
		public void Write (string sequence)
		{
			string segment;
			ulong encoding;

			for (int i = 0; i < sequence.Length; i++) {
				buffer [buffered++] = sequence [i];
				if (buffered == segmentLength) {
					segment = new string (buffer);
					encoding = EncodeSegment (segment);
					buffered = 0;

					WriteBinary (encoding);
				}
			}
		}

		/// <summary>
		/// Writes an encoding result to the file as separate bytes.
		/// </summary>
		/// <param name="encoding">The encoding of a segment</param>
		private void WriteBinary (ulong encoding)
		{
			byte[] bytes = null;
			switch (bits) {
			case 8:
				bytes = new byte[] { (byte)encoding };
				break;

			case 16:
				bytes = BitConverter.GetBytes ((ushort)encoding);
				break;

			case 32:
				bytes = BitConverter.GetBytes ((uint)encoding);
				break;

			case 64:
				bytes = BitConverter.GetBytes (encoding);
				break;
			}

			writer.BaseStream.Write (bytes, 0, bytes.Length);
			writtenSegments++;
		}

		/// <summary>
		/// Encodes what is left in the buffer, writes it to the file, and
		/// closes the reader.
		/// </summary>
		public void End ()
		{
			if (buffered != 0) {
				string segment = new string (buffer, 0, buffered);
				ulong encoding = EncodeSegment (segment);
				buffered = 0;

				WriteBinary (encoding);
			}

			writer.Flush ();
			writer.Close ();
		}

		private ulong EncodeSegment (string segment)
		{
			return GLexEncoding.LexicographicalEncoding (segment, segment.Length);
		}
	}
}

