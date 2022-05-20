import ammonite.ops._

/*
Decorates data files with ids and names.

use,microarray: amm decorateFile.sc geneIdFile=<file:gene id to array design>, geneNameFile=<file:gene id to gene name> source=<file to decorate>
use,rnaseq: amm decorateFile.sc geneNameFile=<file:gene id to gene name> source=<file to decorate>
stdout: decorated file, with amended header.

*/
def fileIntoLines(inFile: Path) = {
	scala.io.Source.fromFile(inFile.toIO)
	.getLines.toStream
	.map(_.split("\t", -1))
}

def fileIntoPairs(inFile: Path) = {
	fileIntoLines(inFile)
	.map{case l => (l.lift(0), l.lift(1))}
	.collect{case (Some(key), Some(value)) => (key, value)}
}

def insertPropertyValueAfterId(map : Map[String, String])(line :Array[String]): Array[String] = {
	Array(line.head, map.getOrElse(line.head, "")) ++ line.tail
}

def tryPrependId(secondColumnToFirstColumn : Map[String, String])(line :Array[String]): Option[Array[String]] = {
	secondColumnToFirstColumn
	.get(line.head)
	.map{Array(_) ++ line}
}

def prependGeneIdLabel(header: Array[String]) = {
	Array("Gene ID") ++ header
}

def addGeneNameLabelAfterIdLabel(header: Array[String]) = {
	Array(header.head, "Gene Name") ++ header.tail
}

def prependGeneId(geneIdFile:Path)(in: => Stream[Array[String]]) = {
	in
	.flatMap(
		tryPrependId(
			fileIntoPairs(geneIdFile)
			.map(_.swap)
			.toMap
		)
	_ )
}

def addGeneNameAfterId(geneNameFile:Path)(in: => Stream[Array[String]]) = {
	in
	.map(
		insertPropertyValueAfterId(
			fileIntoPairs(geneNameFile)
			.toMap
		)
	_ )
}

def transformWhole[T](transformHeader: T => T,transformTail: (=>Stream[T]) => Stream[T])(in: => Stream[T]) = {
	transformHeader(in.head) #:: transformTail(in.tail)
}

def callByName_andThen[X](a: (=> X) => X , b: (=>X) => X)(x: => X) = {
	b(a(x))
}

// https://github.com/scopt/scopt/issues/145
implicit def optionPathRead: scopt.Read[Option[Path]] = scopt.Read.stringRead.map{case "" => None ; case str => Some(Path(str, pwd))}

@main
def decorate(geneIdFile:Option[Path]=None, geneNameFile:Path, source:Path) : Unit = {

	(geneIdFile match {
		case None
			=> transformWhole(addGeneNameLabelAfterIdLabel _ , addGeneNameAfterId(geneNameFile) _ ) _
		case Some(geneIdFile)
			=> transformWhole (
				prependGeneIdLabel _ andThen addGeneNameLabelAfterIdLabel _ ,
				callByName_andThen[Stream[Array[String]]](prependGeneId(geneIdFile) _, addGeneNameAfterId(geneNameFile) _)
			) _
	})(fileIntoLines(source))
	.map(_.mkString("\t"))
	.foreach {
		System.out.println _
	}
}