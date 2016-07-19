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
	/// A reader class for GLex files, currently only supports DNA and only raw A,C,G,T strings.
	/// </summary>
	public class GLexReader
	{
		StreamReader reader;

		// information from the header
		byte segmentLength = 0;
		ulong sequenceLength = 0;
		byte bits;

		// helper variables
		int readSegments = 0;
		int numSegments = 0;

		// buffer-related variables
		char[] buffer;
		long startPosition = 0;
		ulong buffered = 0;
		ulong readLength = 0;

		/// <summary>
		/// Initializes a new instance of the <see cref="GLex.GLexReader"/>.
		/// </summary>
		/// <param name="file">The path to the file which has the GLex-encoded genome</param>
		public GLexReader (string file)
		{
			reader = new StreamReader (file);
			ReadHeader ();
			numSegments = (int)Math.Ceiling(sequenceLength / (double)segmentLength);
		}

		/// <summary>
		/// Reads the header of the GLex file and sets the reader settings accordingly.
		/// </summary>
		private void ReadHeader ()
		{
			bits = (byte)reader.BaseStream.ReadByte ();

			byte[] bytes = new byte[8];
			reader.BaseStream.Read (bytes, 0, bytes.Length);

			sequenceLength = (ulong)BitConverter.ToUInt64 (bytes, 0);
			segmentLength = (byte)reader.BaseStream.ReadByte ();

			// for buffering the stream
			buffer = new char[segmentLength];
			startPosition = reader.BaseStream.Position;
		}

		/// <summary>
		/// Reads and decodes the next segment directly from the file reader.
		/// This function does not take into consideration what is in the buffer,
		/// and does not write anything to the buffer.
		/// </summary>
		/// <returns>The next segment.</returns>
		private string ReadNextSegment ()
		{
			byte[] buffer = new byte[bits / 8];
			ulong encoding = 0;
			string segment;

			if (reader.BaseStream.Read (buffer, 0, buffer.Length) == 0)
				return null;

			switch (bits) {
			case 16:
				encoding = BitConverter.ToUInt16 (buffer, 0);
				break;

			case 32:
				encoding = BitConverter.ToUInt32 (buffer, 0);
				break;

			case 64:
				encoding = BitConverter.ToUInt64 (buffer, 0);
				break;

			default:
				throw new Exception ("Invalid mode");
			}

			readSegments++;

			if (readSegments != numSegments)
				segment = GLexEncoding.LexicographicalDecoding (encoding, segmentLength);
			else {
				segment = GLexEncoding.LexicographicalDecoding (encoding, (int)sequenceLength - (readSegments-1)*segmentLength);
			}


			readLength += (ulong)segment.Length;
			return segment;
		}

		#region BUFFER_FUNCTIONS
		int CopyFromBuffer(char[] destination, int offset, int len)
		{
			if ((ulong)len > buffered)
				throw new ArgumentException("The length of the current buffer is less than the given length");

			for (int i = 0; i < len; i++)
			{
				destination[offset + i] = buffer[i];
				buffered--;
			}

			// shift the buffer to the left
			for (int i = 0; i < buffer.Length - len; i++)
			{
				buffer[i] = buffer[i + len];
			}

			return len;
		}

		int WriteToBuffer(char[] source, int start)
		{
			int initial = (int)buffered;
			for (int i = start; i < source.Length; i++)
			{
				buffer[buffered++] = source[i];
			}
			return (int)buffered-initial;
		}

		int WriteToBuffer(string source, int start)
		{
			int initial = (int)buffered;
			for (int i = start; i < source.Length; i++)
			{
				buffer[buffered++] = source[i];
			}
			return (int)buffered - initial;
		}
		#endregion

		/// <summary>
		/// Reads a segment with the default segmentLength specifed in the file header
		/// </summary>
		/// <returns>The segment.</returns>
		public string ReadSegment ()
		{
			return ReadSegment (segmentLength);
		}

		/// <summary>
		/// Reads a segment of length k
		/// </summary>
		/// <returns>The segment.</returns>
		/// <param name="k">K.</param>
		public string ReadSegment (int k)
		{
			return ReadSegment ((ulong)k);
		}

		/// <summary>
		/// Reads a segment of length k
		/// </summary>
		/// <returns>The segment.</returns>
		/// <param name="k">K.</param>
		public string ReadSegment (ulong k)
		{
			if (k == 0)
				return "";

			if (readLength >= sequenceLength && buffered < 1)
				return null;

			// initialization of variables needed to read the segment
			int length = k <= buffered ?
				(int)Math.Min (k, buffered) :
				(int)Math.Min (k, buffered + sequenceLength - readLength);
			int fromBuffer = Math.Min (length, (int)buffered);
			char[] segment = new char[length];
			int read = 0;
			int segmentIndex = 0; // dictates where the last reading from a segment was, used to copy what is left to the buffer
			string nextSegment = ""; // the fixed-size segments (of size segmentLength from the header)

			// read what is needed from the buffer
			if (fromBuffer > 0)
				read += CopyFromBuffer (segment, read, fromBuffer);

			// keep reading new fixed-size segments until we have a segment of the wanted length
			while (read != length) {
				nextSegment = ReadNextSegment ();
				for (segmentIndex = 0; segmentIndex < nextSegment.Length && read < length; segmentIndex++)
					segment [read++] = nextSegment [segmentIndex];
			}

			// if there is a part of the last fixed-size segment left, write it to the buffer
			if (segmentIndex < nextSegment.Length)
				segmentIndex += WriteToBuffer (nextSegment, segmentIndex);

			return new string (segment);
		}

		/// <summary>
		/// Resets the reader to the position of the data, clears the buffer, and reads the entire genome sequence
		/// </summary>
		/// <returns>The sequence.</returns>
		public string ReadSequence ()
		{
			byte[] buffer = new byte[bits / 8];
			char[] sequence = new char[sequenceLength];
			int index = 0;
			ulong encoding = 0;
			string segment;

			// reset
			reader.BaseStream.Seek (startPosition, SeekOrigin.Begin);
			buffered = 0;
			readLength = 0;
			readSegments = 0;

			while (reader.BaseStream.Read(buffer, 0, buffer.Length) != 0) {
				switch (bits) {
				case 16:
					encoding = BitConverter.ToUInt16 (buffer, 0);
					break;

				case 32:
					encoding = BitConverter.ToUInt32 (buffer, 0);
					break;

				case 64:
					encoding = BitConverter.ToUInt64 (buffer, 0);
					break;

				default:
					throw new Exception ("Invalid mode");
				}
				readSegments++;

				if (readSegments != numSegments)
					segment = GLexEncoding.LexicographicalDecoding (encoding, segmentLength);
				else {
					segment = GLexEncoding.LexicographicalDecoding (encoding, (int)sequenceLength - (readSegments-1)*segmentLength);
				}

				for (int i = 0; i < segment.Length; i++) {
					sequence [index++] = segment [i];
				}

				if (readSegments == numSegments)
					break;
			}

			return new string (sequence);
		}
	}
}

